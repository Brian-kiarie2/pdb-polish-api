FROM mambaorg/micromamba:1.5.8
ARG MAMBA_USER_ID=1000
ARG MAMBA_USER_GID=1000
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:$PATH
WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER app.py .
ENV PORT=7860
EXPOSE 7860
CMD ["sh", "-c", "uvicorn app:app --host 0.0.0.0 --port ${PORT:-7860}"]
