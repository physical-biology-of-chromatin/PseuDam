nextflow src/nf_modules/SAMtools/tests/sort_bams.nf \
  -c src/nf_modules/SAMtools/samtools.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.bam"

nextflow src/nf_modules/SAMtools/tests/index_bams.nf \
  -c src/nf_modules/SAMtools/samtools.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.sort.bam"

nextflow src/nf_modules/SAMtools/tests/split_bams.nf \
  -c src/nf_modules/SAMtools/samtools.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.bam"

nextflow src/nf_modules/SAMtools/tests/filter_bams.nf \
  -c src/nf_modules/SAMtools/samtools.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.bam" \
  --bed "data/tiny_dataset/OLD/2genes.bed"
