./nextflow src/nf_modules/sambamba/sort_bams.nf \
  -c src/nf_modules/sambamba/sort_bams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.bam" \
  -resume

./nextflow src/nf_modules/sambamba/index_bams.nf \
  -c src/nf_modules/sambamba/index_bams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.sort.bam" \
  -resume

./nextflow src/nf_modules/sambamba/split_bams.nf \
  -c src/nf_modules/sambamba/split_bams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.bam" \
  -resume

if [ -x "$(command -v singularity)" ]; then
./nextflow src/nf_modules/sambamba/sort_bams.nf \
  -c src/nf_modules/sambamba/sort_bams.config \
  -profile singularity \
  --bam "data/tiny_dataset/map/tiny_v2.bam" \
  -resume

./nextflow src/nf_modules/sambamba/index_bams.nf \
  -c src/nf_modules/sambamba/index_bams.config \
  -profile singularity \
  --bam "data/tiny_dataset/map/tiny_v2.sort.bam" \
  -resume

./nextflow src/nf_modules/sambamba/split_bams.nf \
  -c src/nf_modules/sambamba/split_bams.config \
  -profile singularity \
  --bam "data/tiny_dataset/map/tiny_v2.bam" \
  -resume
fi
