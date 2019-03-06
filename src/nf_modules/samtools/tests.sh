./nextflow src/nf_modules/samtools/sort_bams.nf \
  -c src/nf_modules/samtools/sort_bams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.bam" \
  -resume

./nextflow src/nf_modules/samtools/index_bams.nf \
  -c src/nf_modules/samtools/index_bams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.sort.bam" \
  -resume

./nextflow src/nf_modules/samtools/split_bams.nf \
  -c src/nf_modules/samtools/split_bams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.bam" \
  -resume

./nextflow src/nf_modules/samtools/filter_bams.nf \
  -c src/nf_modules/samtools/filter_bams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.bam" \
  --bed "data/tiny_dataset/OLD/2genes.bed" \
  -resume

if [ -x "$(command -v singularity)" ]; then
./nextflow src/nf_modules/samtools/sort_bams.nf \
  -c src/nf_modules/samtools/sort_bams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.bam" \
  -resume

./nextflow src/nf_modules/samtools/index_bams.nf \
  -c src/nf_modules/samtools/index_bams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.sort.bam" \
  -resume

./nextflow src/nf_modules/samtools/split_bams.nf \
  -c src/nf_modules/samtools/split_bams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.bam" \
  -resume

./nextflow src/nf_modules/samtools/filter_bams.nf \
  -c src/nf_modules/samtools/filter_bams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.bam" \
  --bed "data/tiny_dataset/OLD/2genes.bed" \
  -resume
fi
