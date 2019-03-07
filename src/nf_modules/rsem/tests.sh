./nextflow src/nf_modules/rsem/indexing.nf \
  -c src/nf_modules/rsem/indexing.config \
  -profile docker \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" \
  --annotation "data/tiny_dataset/annot/tiny.gff" \
  -resume

./nextflow src/nf_modules/rsem/quantification_single.nf \
  -c src/nf_modules/rsem/quantification_single.config \
  -profile docker \
  --index "results/mapping/index/tiny_v2.index*" \
  --fastq "data/tiny_dataset/fastq/tiny*_R1.fastq" \
  -resume

./nextflow src/nf_modules/rsem/quantification_paired.nf \
  -c src/nf_modules/rsem/quantification_paired.config \
  -profile docker \
  --index "results/mapping/index/tiny_v2.index*" \
  --fastq "data/tiny_dataset/fastq/tiny*_R{1,2}.fastq" \
  -resume


if [ -x "$(command -v singularity)" ]; then
./nextflow src/nf_modules/rsem/indexing.nf \
  -c src/nf_modules/rsem/indexing.config \
  -profile docker \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" \
  --annotation "data/tiny_dataset/annot/tiny.gff" \
  -resume

./nextflow src/nf_modules/rsem/quantification_single.nf \
  -c src/nf_modules/rsem/quantification_single.config \
  -profile docker \
  --index "results/mapping/index/tiny_v2.index*" \
  --fastq "data/tiny_dataset/fastq/tiny*_R1.fastq" \
  -resume

./nextflow src/nf_modules/rsem/quantification_paired.nf \
  -c src/nf_modules/rsem/quantification_paired.config \
  -profile docker \
  --index "results/mapping/index/tiny_v2.index*" \
  --fastq "data/tiny_dataset/fastq/tiny*_R{1,2}.fastq" \
  -resume
fi
