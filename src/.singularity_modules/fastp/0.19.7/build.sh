#/bin/sh
sudo singularity build --force bin/fastp:0.19.7.img src/singularity_modules/fastp/0.19.7/fastp.def && \
singularity sign bin/fastp:0.19.7.img
