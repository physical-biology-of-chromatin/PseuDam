#/bin/sh
sudo singularity build --force bin/deeptools:3.0.2.img src/singularity_modules/deeptools/3.0.2/deeptools.def && \
singularity sign bin/deeptools:3.0.2.img