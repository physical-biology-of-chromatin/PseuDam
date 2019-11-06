#/bin/sh
sudo singularity build --force bin/music:6613c53.img src/singularity_modules/music/6613c53/music.def && \
singularity sign bin/music:6613c53.img
