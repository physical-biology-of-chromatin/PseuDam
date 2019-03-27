#/bin/sh
sudo singularity build --force bin/hisat2:2.0.0.sif src/singularity_modules/hisat2/2.0.0/hisat2.def && \
singularity sign bin/hisat2:2.0.0.sif
