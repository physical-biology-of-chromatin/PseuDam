#/bin/sh
sudo singularity build --force bin/subread:1.32.2.sif docker://quay.io/biocontainers/subread:1.6.4--h84994c4_1 && \
singularity sign bin/subread:1.32.2.sif
