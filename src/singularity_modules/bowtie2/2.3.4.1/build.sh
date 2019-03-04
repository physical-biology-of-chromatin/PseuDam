#/bin/sh
sudo singularity build --force bin/bowtie2:2.3.4.1.sif src/singularity_modules/Bowtie2/2.3.4.1/Bowtie2.def
singularity sign bin/bowtie2:2.3.4.1.sif
