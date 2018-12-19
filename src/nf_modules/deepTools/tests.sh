#!/bin/sh

./nextflow src/nf_modules/deepTools/bam_to_bigwig.nf -c src/nf_modules/deepTools/bam_to_bigwig.config -profile docker --bam "data/tiny_dataset/map/tiny_v2.sort.bam"

./nextflow src/nf_modules/deepTools/compute_matrix.nf -c src/nf_modules/deepTools/compute_matrix.config -profile docker --bam "results/mapping/bigwig/*.bw"

