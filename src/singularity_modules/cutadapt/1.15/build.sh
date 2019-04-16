#/bin/sh
sudo singularity build --force bin/cutadapt:1.15.img src/singularity_modules/cutadapt/1.15/cutadapt.def && \
singularity sign bin/cutadapt:1.15.img
