#/bin/sh
sudo singularity build --force bin/htseq:0.8.0.sif src/singularity_modules/htseq/0.8.0/htseq.def && \
singularity sign bin/htseq:0.8.0.sif
