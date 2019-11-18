#/bin/sh
sudo singularity build --force bin/fastqc:0.11.5.img src/singularity_modules/fastqc/0.11.5/fastqc.def && \
singularity sign bin/fastqc:0.11.5.img
