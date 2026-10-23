"""pdb-polish API: missing atoms, polar/all hydrogens (PDBFixer) and a protein-only PDBQT (OpenBabel)."""
import os
import shutil
import subprocess
import tempfile
import threading
from io import StringIO

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from openmm import Context, CustomExternalForce, LocalEnergyMinimizer, Platform, VerletIntegrator, unit
from openmm import app as omapp
from pdbfixer import PDBFixer

MAX_BYTES = 15 * 1024 * 1024
slots = threading.BoundedSemaphore(2)  # at most 2 jobs at once on a small free server
AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HID HIE HIP CYX ASH GLH LYN MSE SEC PYL".split())

app = FastAPI(title="pdb-polish API")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])


@app.get("/health")
def health():
    return {"ok": True, "obabel": shutil.which("obabel") is not None}


def relax_hydrogens(fixer, log):
    """Short minimization with all heavy atoms restrained, so only hydrogens move."""
    ff = omapp.ForceField("amber14-all.xml", "amber14/tip3p.xml")
    system = ff.createSystem(fixer.topology, nonbondedMethod=omapp.NoCutoff)
    force = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter("k", 50000.0)
    for name in ("x0", "y0", "z0"):
        force.addPerParticleParameter(name)
    for i, (atom, pos) in enumerate(zip(fixer.topology.atoms(), fixer.positions)):
        if atom.element is not None and atom.element.symbol != "H":
            v = pos.value_in_unit(unit.nanometer)
            force.addParticle(i, [v[0], v[1], v[2]])
    system.addForce(force)
    ctx = Context(system, VerletIntegrator(0.001 * unit.picoseconds), Platform.getPlatformByName("CPU"))
    ctx.setPositions(fixer.positions)
    LocalEnergyMinimizer.minimize(ctx, 10.0, 300)
    fixer.positions = ctx.getState(getPositions=True).getPositions()
    log.append("hydrogens relaxed (heavy atoms restrained)")


def strip_nonpolar_h(topology, positions):
    """Delete hydrogens bonded to carbon; keep those on N, O, S."""
    mod = omapp.Modeller(topology, positions)
    drop = []
    for a1, a2 in mod.topology.bonds():
        for h, other in ((a1, a2), (a2, a1)):
            if h.element is not None and h.element.symbol == "H" and other.element is not None and other.element.symbol == "C":
                drop.append(h)
    mod.delete(drop)
    return mod.topology, mod.positions, len(drop)


def protein_only(topology, positions):
    mod = omapp.Modeller(topology, positions)
    extra = [r for r in mod.topology.residues() if r.name not in AA]
    mod.delete(extra)
    return mod.topology, mod.positions, len(extra)


def count_h(topology):
    return sum(1 for a in topology.atoms() if a.element is not None and a.element.symbol == "H")


@app.post("/prepare")
def prepare(
    file: UploadFile = File(...),
    ph: float = Form(7.4),
    add_atoms: bool = Form(True),
    add_h: bool = Form(True),
    h_mode: str = Form("polar"),
    replace_nonstd: bool = Form(False),
    minimize: bool = Form(False),
    pdbqt: bool = Form(True),
):
    if not 0 <= ph <= 14:
        raise HTTPException(422, "pH must be between 0 and 14.")
    if h_mode not in ("polar", "all"):
        raise HTTPException(422, "h_mode must be 'polar' or 'all'.")
    data = file.file.read(MAX_BYTES + 1)
    if len(data) > MAX_BYTES:
        raise HTTPException(413, "File is larger than 15 MB.")
    if b"ATOM" not in data and b"HETATM" not in data:
        raise HTTPException(422, "No ATOM/HETATM records found.")
    if not slots.acquire(blocking=False):
        raise HTTPException(503, "Server is busy. Try again in a minute.")
    try:
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "in.pdb")
            with open(src, "wb") as fh:
                fh.write(data)
            fixer = PDBFixer(filename=src)
            log = []
            fixer.findMissingResidues()
            fixer.missingResidues = {}  # gaps are detected by the page but not modelled
            fixer.findNonstandardResidues()
            if replace_nonstd:
                log.append(f"{len(fixer.nonstandardResidues)} nonstandard residues replaced")
                fixer.replaceNonstandardResidues()
            fixer.findMissingAtoms()
            if add_atoms:
                n = sum(len(v) for v in fixer.missingAtoms.values()) + sum(len(v) for v in fixer.missingTerminals.values())
                fixer.addMissingAtoms()
                log.append(f"{n} missing heavy atoms added")
            topology, positions = fixer.topology, fixer.positions
            if add_h:
                fixer.addMissingHydrogens(ph)
                log.append(f"hydrogens added for pH {ph}")
                if minimize:  # needs all hydrogens, so it runs before any polar-only filtering
                    try:
                        relax_hydrogens(fixer, log)
                    except Exception as e:  # e.g. ligand without a force-field template
                        log.append(f"minimization skipped: {str(e)[:120]}")
                topology, positions = fixer.topology, fixer.positions
                if h_mode == "polar":
                    topology, positions, n = strip_nonpolar_h(topology, positions)
                    log.append(f"polar hydrogens kept; {n} carbon-bound hydrogens removed")
                log.append(f"{count_h(topology)} hydrogens in output")
            buf = StringIO()
            omapp.PDBFile.writeFile(topology, positions, buf, keepIds=True)
            pdb_text = buf.getvalue()
            qt_text = None
            if pdbqt:
                obabel = shutil.which("obabel")
                if obabel is None:
                    log.append("PDBQT unavailable: obabel not found on the server")
                else:
                    ptop, ppos, removed = protein_only(topology, positions)
                    if ptop.getNumAtoms() == 0:
                        log.append("PDBQT skipped: no protein atoms left")
                    else:
                        out_pdb, out_qt = os.path.join(d, "rec.pdb"), os.path.join(d, "rec.pdbqt")
                        with open(out_pdb, "w") as fh:
                            omapp.PDBFile.writeFile(ptop, ppos, fh, keepIds=True)
                        r = subprocess.run([obabel, out_pdb, "-O", out_qt, "-xr"], capture_output=True, text=True, timeout=120)
                        if os.path.exists(out_qt) and os.path.getsize(out_qt) > 0:
                            with open(out_qt) as fh:
                                qt_text = fh.read()
                            log.append("PDBQT written for the protein only (rigid receptor, OpenBabel Gasteiger charges)")
                            if removed:
                                log.append(f"{removed} non-protein residues (waters, ligands, ions) left out of the PDBQT")
                            if not add_h:
                                log.append("warning: no hydrogens were added, so the PDBQT has no polar H")
                        else:
                            log.append("PDBQT failed: " + (r.stderr or "unknown error")[:120])
            return {"pdb": pdb_text, "pdbqt": qt_text, "log": log}
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(422, f"Could not prepare this structure: {str(e)[:200]}. If heteroatoms are present, remove them in the Heteroatoms tab first.")
    finally:
        slots.release()
