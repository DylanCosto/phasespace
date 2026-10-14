# PhaseSpace container: the CLI and its only dependency (pysam), no network use at runtime.
# Build: docker build -t phasespace .   Run: docker run --rm -v "$PWD:/data" -w /data phasespace phasespace --help
FROM python:3.12-slim

LABEL org.opencontainers.image.title="PhaseSpace" \
      org.opencontainers.image.description="Explore phase uncertainty and joint coding consequences in VCF/gVCF data" \
      org.opencontainers.image.source="https://github.com/DylanCosto/phasespace" \
      org.opencontainers.image.licenses="MIT"

# procps provides ps, which Nextflow uses to collect task metrics.
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/phasespace
COPY pyproject.toml README.md MANIFEST.in LICENSE ./
COPY src ./src
RUN pip install --no-cache-dir . && phasespace --version

# No ENTRYPOINT: workflow managers run shell commands inside the image.
CMD ["phasespace", "--help"]
