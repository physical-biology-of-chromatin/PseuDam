nextflow src/nf_modules/FastQC/tests/fastqc_paired.nf \
  -c src/nf_modules/FastQC/fastqc.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq"

nextflow src/nf_modules/FastQC/tests/fastqc_single.nf \
  -c src/nf_modules/FastQC/fastqc.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_S.fastq"
