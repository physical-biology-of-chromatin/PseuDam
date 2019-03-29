./nextflow src/nf_modules/hisat2/indexing.nf \
  -c src/nf_modules/hisat2/indexing.config \
  -profile docker \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" \
  -resume

./nextflow src/nf_modules/hisat2/mapping_paired.nf \
  -c src/nf_modules/hisat2/mapping_paired.config \
  -profile docker \
  --index "results/mapping/index/tiny_v2.index*" \
  --fastq "data/tiny_dataset/fastq/tiny*_R{1,2}.fastq" \
  -resume

./nextflow src/nf_modules/hisat2/mapping_single.nf \
  -c src/nf_modules/hisat2/mapping_single.config \
  -profile docker \
  --index "results/mapping/index/tiny_v2.index*" \
  --fastq "data/tiny_dataset/fastq/tiny*_S.fastq" \
  -resume

if [ -x "$(command -v singularity)" ]; then
./nextflow src/nf_modules/hisat2/indexing.nf \
  -c src/nf_modules/hisat2/indexing.config \
  -profile singularity \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" \
  -resume

./nextflow src/nf_modules/hisat2/mapping_paired.nf \
  -c src/nf_modules/hisat2/mapping_paired.config \
  -profile singularity \
  --index "results/mapping/index/tiny_v2.index*" \
  --fastq "data/tiny_dataset/fastq/tiny*_R{1,2}.fastq"

./nextflow src/nf_modules/hisat2/mapping_single.nf \
  -c src/nf_modules/hisat2/mapping_single.config \
  -profile singularity \
  --index "results/mapping/index/tiny_v2.index*" \
  --fastq "data/tiny_dataset/fastq/tiny*_S.fastq"
fi
