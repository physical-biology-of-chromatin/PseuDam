#/bin/sh
sudo singularity build --force bin/htseq:0.11.2.sif src/singularity_modules/htseq/0.11.2/htseq.def && \
singularity sign bin/htseq:0.11.2.sif
