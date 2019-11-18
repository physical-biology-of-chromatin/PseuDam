#/bin/sh
sudo singularity build --force bin/multiqc:1.7.img src/singularity_modules/multiqc/1.7/multiqc.def && \
singularity sign bin/multiqc:1.7.img
