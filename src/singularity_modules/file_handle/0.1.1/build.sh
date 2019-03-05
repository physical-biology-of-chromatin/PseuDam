#/bin/sh
sudo singularity build --force bin/file_handle:0.1.1.sif src/singularity_modules/file_handle/0.1.1/file_handle.def && \
singularity sign bin/file_handle:0.1.1.sif
