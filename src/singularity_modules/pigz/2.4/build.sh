#/bin/sh
sudo singularity build --force bin/pigz:2.4.sif src/singularity_modules/pigz/2.4/pigz.def
singularity sign bin/pigz:2.4.sif
