#/bin/sh
sudo singularity build --force bin/bedtools:2.25.0.sif src/singularity_modules/BEDtools/2.25.0/BEDtools.def
singularity sign bin/bedtools:2.25.0.sif
