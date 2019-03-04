#/bin/sh
sudo singularity build --force bin/bwa:0.7.17.sif src/singularity_modules/bwa/0.7.17/bwa.def
singularity sign bin/bwa:0.7.17.sif
