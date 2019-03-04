#/bin/sh
sudo singularity build --force bin/sratoolkit:2.8.2.sif src/singularity_modules/sratoolkit/2.8.2/sratoolkit.def
singularity sign bin/sratoolkit:2.8.2.sif
