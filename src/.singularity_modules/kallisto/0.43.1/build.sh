#/bin/sh
sudo singularity build --force bin/kallisto:0.43.1.img src/singularity_modules/kallisto/0.43.1/kallisto.def && \
singularity sign bin/kallisto:0.43.1.img