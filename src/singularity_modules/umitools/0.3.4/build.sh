#/bin/sh
sudo singularity build --force bin/umitools:0.3.4.sif docker://quay.io/biocontainers/umitools:0.3.4--py37_1 && \
singularity sign bin/umitools:0.3.4.sif
