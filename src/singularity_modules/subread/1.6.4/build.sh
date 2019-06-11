#/bin/sh
sudo singularity build --force bin/subread:1.6.4.img src/singularity_modules/subread/1.6.4/subread.def && \
singularity sign bin/subread:1.6.4.img
