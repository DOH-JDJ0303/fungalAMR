FROM mambaorg/micromamba:1.5.8
USER root
ENV DEBIAN_FRONTEND=noninteractive

# Install tools into the base env
RUN micromamba install -y -n base -c conda-forge -c bioconda --channel-priority strict \
      python=3.11 freebayes=1.3.6 bcftools=1.20 bwa=0.7.18 samtools=1.20 htslib=1.20 seqtk=1.4 \
  && micromamba clean -a -y

# App
WORKDIR /app
COPY --chown=mambauser:mambauser bin/ /app/
RUN chmod -R +x /app

# Make the env “active” by default
ENV PATH="/opt/conda/bin:$PATH" \
    HOME="/home/mambauser" \
    XDG_CACHE_HOME="/home/mambauser/.cache"
RUN mkdir -p /home/mambauser && chown -R mambauser:mambauser /home/mambauser

# Default user & workdir
USER mambauser
WORKDIR /data

ENTRYPOINT ["python", "/app/fungalAMR"]
