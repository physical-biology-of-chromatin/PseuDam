#/bin/sh
sudo singularity build --force bin/umi_tools:1.0.0.sif src/singularity_modules/umi_tools/1.0.0/umi_tools.def && \
singularity sign bin/umi_tools:1.0.0.sif
