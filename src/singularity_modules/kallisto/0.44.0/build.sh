#/bin/sh
sudo singularity build --force bin/kallisto:0.44.0.sif src/singularity_modules/kallisto/0.44.0/kallisto.def && \
singularity sign bin/kallisto:0.44.0.sif
