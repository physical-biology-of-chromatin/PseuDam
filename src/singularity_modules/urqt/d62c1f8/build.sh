#/bin/sh
sudo singularity build --force bin/urqt:d62c1f8.img src/singularity_modules/urqt/d62c1f8/urqt.def && \
singularity sign bin/urqt:d62c1f8.img
