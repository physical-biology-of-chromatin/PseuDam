#/bin/sh
sudo singularity build --force bin/tophat:2.1.1.sif src/singularity_modules/tophat/2.1.1/tophat.def && \
singularity sign bin/tophat:2.1.1.sif
