#/bin/sh
sudo singularity build --force bin/cutadapt:2.1.img src/singularity_modules/cutadapt/2.1/cutadapt.def && \
singularity sign bin/cutadapt:2.1.img
