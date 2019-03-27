./nextflow src/nf_modules/bwa/indexing.nf \
  -c src/nf_modules/bwa/indexing.config \
  -profile docker \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" \
  -resume

# ./nextflow src/nf_modules/bwa/mapping_single.nf \
#   -c src/nf_modules/bwa/mapping_single.config \
#   -profile docker \
#   --index "results/mapping/index/tiny_v2.index" \
#   --fastq "data/tiny_dataset/fastq/tiny*_S.fastq"

./nextflow src/nf_modules/bwa/mapping_paired.nf \
  -c src/nf_modules/bwa/mapping_paired.config \
  -profile docker \
  --index "results/mapping/index/tiny_v2.*" \
  --fastq "data/tiny_dataset/fastq/tiny*_R{1,2}.fastq" \
  -resume


if [ -x "$(command -v singularity)" ]; then
./nextflow src/nf_modules/bwa/indexing.nf \
  -c src/nf_modules/bwa/indexing.config \
  -profile singularity \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" \
  -resume


# ./nextflow src/nf_modules/bwa/mapping_single.nf \
#   -c src/nf_modules/bwa/mapping_single.config \
#   -profile singularity \
#   --index "results/mapping/index/tiny_v2.index" \
#   --fastq "data/tiny_dataset/fastq/tiny*_S.fastq"

./nextflow src/nf_modules/bwa/mapping_paired.nf \
  -c src/nf_modules/bwa/mapping_paired.config \
  -profile singularity \
  --index "results/mapping/index/tiny_v2.*" \
  --fastq "data/tiny_dataset/fastq/tiny*_R{1,2}.fastq" \
  -resume

fi
