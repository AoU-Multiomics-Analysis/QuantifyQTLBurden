## Use micromamba as the base image
FROM mambaorg/micromamba:1.5.3
#FROM ghcr.io/prefix-dev/pixi:latest

## Set up environment
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN micromamba install -y -n base -c conda-forge -c bioconda \
    bioconda::htslib \
    bioconda::bcftools
CMD ["bash"]
