#/bin/sh
sudo singularity build --force bin/rsem:1.3.0.sif src/singularity_modules/rsem/1.3.0/rsem.def && \
singularity sign bin/rsem:1.3.0.sif
