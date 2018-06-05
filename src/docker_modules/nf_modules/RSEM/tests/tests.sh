nextflow src/nf_modules/RSEM/tests/index.nf \
  -c src/nf_modules/RSEM/rsem.config \
  -profile docker \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" \
  --annotation "data/tiny_dataset/annot/tiny.gff"

nextflow src/nf_modules/RSEM/tests/quantification_single.nf \
  -c src/nf_modules/RSEM/rsem.config \
  -profile docker \
  --index "results/mapping/index/tiny_v2.index*" \
  --fastq "data/tiny_dataset/fastq/tiny*_S.fastq"

nextflow src/nf_modules/RSEM/tests/quantification_paired.nf \
  -c src/nf_modules/RSEM/rsem.config \
  -profile docker \
  --index "results/mapping/index/tiny_v2.index*" \
  --fastq "data/tiny_dataset/fastq/tiny*_R{1,2}.fastq"

