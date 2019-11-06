#/bin/sh
sudo singularity build --force bin/gatk:4.0.8.1.img src/singularity_modules/gatk/4.0.8.1/gatk.def && \
singularity sign bin/gatk:4.0.8.1.img
