#/bin/sh
sudo singularity build --force bin/bowtie:1.2.2.sif src/singularity_modules/bowtie/1.2.2/bowtie.def && \
singularity sign bin/bowtie:1.2.2.sif
