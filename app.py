"""ProteinPolish API: adds missing atoms/hydrogens (PDBFixer) and makes PDBQT (OpenBabel)."""
import os
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

app = FastAPI(title="ProteinPolish API")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])


@app.get("/health")
def health():
    return {"ok": True}


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


@app.post("/prepare")
def prepare(
    file: UploadFile = File(...),
    ph: float = Form(7.4),
    add_atoms: bool = Form(True),
    add_h: bool = Form(True),
    replace_nonstd: bool = Form(False),
    minimize: bool = Form(False),
    pdbqt: bool = Form(True),
):
    if not 0 <= ph <= 14:
        raise HTTPException(422, "pH must be between 0 and 14.")
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
            fixer.missingResidues = {}  # loop modelling is not offered (needs SEQRES)
            fixer.findNonstandardResidues()
            if replace_nonstd:
                log.append(f"{len(fixer.nonstandardResidues)} nonstandard residues replaced")
                fixer.replaceNonstandardResidues()
            fixer.findMissingAtoms()
            if add_atoms:
                n = sum(len(v) for v in fixer.missingAtoms.values()) + sum(len(v) for v in fixer.missingTerminals.values())
                fixer.addMissingAtoms()
                log.append(f"{n} missing heavy atoms added")
            if add_h:
                fixer.addMissingHydrogens(ph)
                log.append(f"hydrogens added for pH {ph}")
                if minimize:
                    try:
                        relax_hydrogens(fixer, log)
                    except Exception as e:  # e.g. ligand without a force-field template
                        log.append(f"minimization skipped: {str(e)[:120]}")
            buf = StringIO()
            omapp.PDBFile.writeFile(fixer.topology, fixer.positions, buf, keepIds=True)
            pdb_text = buf.getvalue()
            qt_text = None
            if pdbqt:
                out_pdb, out_qt = os.path.join(d, "out.pdb"), os.path.join(d, "out.pdbqt")
                with open(out_pdb, "w") as fh:
                    fh.write(pdb_text)
                r = subprocess.run(["obabel", out_pdb, "-O", out_qt, "-xr"], capture_output=True, text=True, timeout=120)
                if os.path.exists(out_qt) and os.path.getsize(out_qt) > 0:
                    with open(out_qt) as fh:
                        qt_text = fh.read()
                    log.append("PDBQT written (rigid receptor, Gasteiger charges from OpenBabel)")
                    if not add_h:
                        log.append("warning: no hydrogens were added, so PDBQT has no polar H")
                else:
                    log.append("PDBQT failed: " + (r.stderr or "unknown error")[:120])
            return {"pdb": pdb_text, "pdbqt": qt_text, "log": log}
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(422, f"Could not prepare this structure: {str(e)[:200]}. If ligands are present, try Heterogens: Water only or Remove all.")
    finally:
        slots.release()
