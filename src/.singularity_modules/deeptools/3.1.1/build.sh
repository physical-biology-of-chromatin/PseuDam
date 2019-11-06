#/bin/sh
sudo singularity build --force bin/deeptools:3.1.1.img src/singularity_modules/deeptools/3.1.1/deeptools.def && \
singularity sign bin/deeptools:3.1.1.img
