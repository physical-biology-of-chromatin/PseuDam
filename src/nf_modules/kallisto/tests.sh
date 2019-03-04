./nextflow src/nf_modules/kallisto/indexing.nf \
  -c src/nf_modules/kallisto/indexing.config \
  -profile docker \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta"

./nextflow src/nf_modules/kallisto/mapping_single.nf \
  -c src/nf_modules/kallisto/mapping_single.config \
  -profile docker \
  --index "results/mapping/index/tiny_v2.index" \
  --fastq "data/tiny_dataset/fastq/tiny*_S.fastq"

./nextflow src/nf_modules/kallisto/mapping_paired.nf \
  -c src/nf_modules/kallisto/mapping_paired.config \
  -profile docker \
  --index "results/mapping/index/tiny_v2.index" \
  --fastq "data/tiny_dataset/fastq/tiny*_R{1,2}.fastq"

