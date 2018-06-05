nextflow src/nf_modules/Kallisto/tests/index.nf \
  -c src/nf_modules/Kallisto/kallisto.config \
  -profile docker \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta"

nextflow src/nf_modules/Kallisto/tests/mapping_single.nf \
  -c src/nf_modules/Kallisto/kallisto.config \
  -profile docker \
  --index "results/mapping/index/tiny_v2.index" \
  --fastq "data/tiny_dataset/fastq/tiny*_S.fastq"

nextflow src/nf_modules/Kallisto/tests/mapping_paired.nf \
  -c src/nf_modules/Kallisto/kallisto.config \
  -profile docker \
  --index "results/mapping/index/tiny_v2.index" \
  --fastq "data/tiny_dataset/fastq/tiny*_R{1,2}.fastq"

