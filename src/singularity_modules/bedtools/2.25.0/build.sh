#/bin/sh
sudo singularity build --force bin/bedtools:2.25.0.sif src/singularity_modules/bedtools/2.25.0/bedtools.def && \
singularity sign bin/bedtools:2.25.0.sif
