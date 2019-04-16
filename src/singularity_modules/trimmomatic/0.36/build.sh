#/bin/sh
sudo singularity build --force bin/trimmomatic:0.36.img src/singularity_modules/trimmomatic/0.36/trimmomatic.def && \
singularity sign bin/trimmomatic:0.36.img
