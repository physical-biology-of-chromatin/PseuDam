#/bin/sh
sudo singularity build --force bin/samblaster:0.1.24.img src/singularity_modules/samblaster/0.1.24/samblaster.def && \
singularity sign bin/samblaster:0.1.24.img