./nextflow src/nf_modules/bowtie2/indexing.nf \
  -c src/nf_modules/bowtie2/indexing.config \
  -profile docker \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" \
  -resume

./nextflow src/nf_modules/bowtie2/mapping_single.nf \
  -c src/nf_modules/bowtie2/mapping_single.config \
  -profile docker \
  --index "data/tiny_dataset/fasta/*.bt2" \
  --fastq "data/tiny_dataset/fastq/tiny*_S.fastq" \
  -resume

./nextflow src/nf_modules/bowtie2/mapping_paired.nf \
  -c src/nf_modules/bowtie2/mapping_paired.config \
  -profile docker \
  --index "data/tiny_dataset/fasta/*.bt2" \
  --fastq "data/tiny_dataset/fastq/tiny*_R{1,2}.fastq" \
  -resume

if [ -x "$(command -v singularity)" ]; then
./nextflow src/nf_modules/bowtie2/indexing.nf \
  -c src/nf_modules/bowtie2/indexing.config \
  -profile singularity \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" \
  -resume

./nextflow src/nf_modules/bowtie2/mapping_single.nf \
  -c src/nf_modules/bowtie2/mapping_single.config \
  -profile singularity \
  --index "data/tiny_dataset/fasta/*.bt2" \
  --fastq "data/tiny_dataset/fastq/tiny*_S.fastq" \
  -resume

./nextflow src/nf_modules/bowtie2/mapping_paired.nf \
  -c src/nf_modules/bowtie2/mapping_paired.config \
  -profile singularity \
  --index "data/tiny_dataset/fasta/*.bt2" \
  --fastq "data/tiny_dataset/fastq/tiny*_R{1,2}.fastq" \
  -resume
fi
