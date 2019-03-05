#/bin/sh
sudo singularity build --force bin/picard:2.18.11.sif src/singularity_modules/picard/2.18.11/picard.def && \
singularity sign bin/picard:2.18.11.sif
