#/bin/sh
sudo singularity build --force bin/sambamba:0.6.7.img src/singularity_modules/sambamba/0.6.7/sambamba.def && \
singularity sign bin/sambamba:0.6.7.img
