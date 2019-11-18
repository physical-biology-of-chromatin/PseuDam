#/bin/sh
sudo singularity build --force bin/samtools:1.7.img src/singularity_modules/samtools/1.7/samtools.def && \
singularity sign bin/samtools:1.7.img
