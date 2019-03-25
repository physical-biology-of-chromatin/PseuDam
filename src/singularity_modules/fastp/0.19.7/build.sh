#/bin/sh
sudo singularity build --force bin/fastp:0.19.7.sif docker://quay.io/biocontainers/fastp:0.19.7--hdbcaa40_0 && \
singularity sign bin/fastp:0.19.7.sif
