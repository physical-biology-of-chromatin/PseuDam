#/bin/sh
sudo singularity build --force bin/bioawk:1.0.sif src/singularity_modules/bioawk/1.0/bioawk.def && \
singularity sign bin/bioawk:1.0.sif
