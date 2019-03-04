#/bin/sh
sudo singularity build --force bin/bcftools:1.7.sif src/singularity_modules/bcftools/1.7/bcftools.def
singularity sign bin/bcftools:1.7.sif
