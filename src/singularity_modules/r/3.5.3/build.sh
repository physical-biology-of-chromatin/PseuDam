#/bin/sh
sudo singularity build --force bin/r:3.5.3.img src/singularity_modules/r/3.5.3/r.def && \
singularity sign bin/r:3.5.3.img
