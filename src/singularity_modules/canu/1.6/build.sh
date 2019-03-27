#/bin/sh
sudo singularity build --force bin/canu:1.6.sif src/singularity_modules/canu/1.6/canu.def && \
singularity sign bin/canu:1.6.sif
