#/bin/sh
sudo singularity build --force bin/ucsc:375.img src/singularity_modules/ucsc/375/ucsc.def && \
singularity sign bin/ucsc:375.img
