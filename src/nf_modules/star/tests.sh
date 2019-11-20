./nextflow src/nf_modules/star/indexing.nf \
  -c src/nf_modules/star/indexing.config \
  -profile docker \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" \
  --annotation "data/tiny_dataset/annot/tiny.gtf" \
  -resume

./nextflow src/nf_modules/star/mapping_single.nf \
  -c src/nf_modules/star/mapping_single.config \
  -profile docker \
  --index "results/mapping/index/*" \
  --fastq "data/tiny_dataset/fastq/tiny*_S.fastq" \
  -resume

./nextflow src/nf_modules/star/mapping_paired.nf \
  -c src/nf_modules/star/mapping_paired.config \
  -profile docker \
  --index "results/mapping/index/*" \
  --fastq "data/tiny_dataset/fastq/tiny*_R{1,2}.fastq" \
  -resume

if [ -x "$(command -v singularity)" ]; then
./nextflow src/nf_modules/star/indexing.nf \
  -c src/nf_modules/star/indexing.config \
  -profile singularity \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" \
  --annotation "data/tiny_dataset/annot/tiny.gtf" \
  -resume

./nextflow src/nf_modules/star/mapping_single.nf \
  -c src/nf_modules/star/mapping_single.config \
  -profile singularity \
  --index "results/mapping/index/*" \
  --fastq "data/tiny_dataset/fastq/tiny*_S.fastq" \
  -resume

./nextflow src/nf_modules/star/mapping_paired.nf \
  -c src/nf_modules/star/mapping_paired.config \
  -profile singularity \
  --index "results/mapping/index/*" \
  --fastq "data/tiny_dataset/fastq/tiny*_R{1,2}.fastq" \
  -resume
fi
