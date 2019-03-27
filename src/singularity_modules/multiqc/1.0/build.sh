#/bin/sh
sudo singularity build --force bin/multiqc:1.0.sif src/singularity_modules/multiqc/1.0/multiqc.def && \
singularity sign bin/multiqc:1.0.sif
