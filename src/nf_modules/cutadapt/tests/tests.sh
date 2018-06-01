nextflow src/nf_modules/cutadapt/tests/adaptor_removal_paired.nf \
  -c src/nf_modules/cutadapt/cutadapt.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq"

nextflow src/nf_modules/cutadapt/tests/adaptor_removal_single.nf \
  -c src/nf_modules/cutadapt/cutadapt.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq"

nextflow src/nf_modules/cutadapt/tests/trimming_paired.nf \
  -c src/nf_modules/cutadapt/cutadapt.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq"

nextflow src/nf_modules/cutadapt/tests/trimming_single.nf \
  -c src/nf_modules/cutadapt/cutadapt.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq"
