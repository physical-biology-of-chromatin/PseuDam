#/bin/sh
sudo singularity build --force bin/cutadapt:1.14.sif src/singularity_modules/cutadapt/1.14/cutadapt.def
singularity sign bin/cutadapt:1.14.sif
