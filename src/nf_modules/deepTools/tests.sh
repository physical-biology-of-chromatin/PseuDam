#!/bin/sh

cp data/tiny_dataset/map/tiny_v2.sort.bam data/tiny_dataset/map/tiny_v2_bis.sort.bam

./nextflow src/nf_modules/deepTools/bam_to_bigwig.nf -c src/nf_modules/deepTools/bam_to_bigwig.config -profile docker --bam "data/tiny_dataset/map/tiny_v2*.sort.bam"

./nextflow src/nf_modules/deepTools/compute_matrix.nf -c src/nf_modules/deepTools/compute_matrix.config -profile docker --bw "results/mapping/bigwig/*.bw" --bed "data/tiny_dataset/annot/tiny.bed"

./nextflow src/nf_modules/deepTools/plot_profile.nf -c src/nf_modules/deepTools/plot_profile.config -profile docker --matrix "results/mapping/region_matrix/*.mat.gz" --title "plot title"
