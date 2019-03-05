#/bin/sh
sudo singularity build --force bin/trimmomatic:0.36.sif src/singularity_modules/trimmomatic/0.36/trimmomatic.def
singularity sign bin/trimmomatic:0.36.sif
