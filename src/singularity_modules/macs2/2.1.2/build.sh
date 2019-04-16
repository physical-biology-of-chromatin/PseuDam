#/bin/sh
sudo singularity build --force bin/macs2:2.1.2.img docker://quay.io/biocontainers/macs2:2.1.2--py27r351h14c3975_1 && \
singularity sign bin/macs2:2.1.2.img
