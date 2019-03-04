#/bin/sh
sudo singularity build --force bin/salmon:0.8.2.sif src/singularity_modules/salmon/0.8.2/salmon.def
singularity sign bin/salmon:0.8.2.sif
