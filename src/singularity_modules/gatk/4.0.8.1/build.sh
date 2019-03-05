#/bin/sh
sudo singularity build --force bin/gatk:4.0.8.1.sif src/singularity_modules/gatk/4.0.8.1/gatk.def && \
singularity sign bin/gatk:4.0.8.1.sif
