#/bin/sh
sudo singularity build --force bin/sratoolkit:2.8.2.img src/singularity_modules/sratoolkit/2.8.2/sratoolkit.def && \
singularity sign bin/sratoolkit:2.8.2.img
