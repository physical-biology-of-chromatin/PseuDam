#/bin/sh
sudo singularity build --force bin/bowtie:1.2.2.sif src/singularity_modules/Bowtie/1.2.2/Bowtie.def
singularity sign bin/bowtie:1.2.2.sif
