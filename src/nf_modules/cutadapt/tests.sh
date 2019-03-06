./nextflow src/nf_modules/cutadapt/adaptor_removal_paired.nf \
  -c src/nf_modules/cutadapt/adaptor_removal_paired.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq" \
  -resume

./nextflow src/nf_modules/cutadapt/adaptor_removal_single.nf \
  -c src/nf_modules/cutadapt/adaptor_removal_single.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny*_S.fastq" \
  -resume

./nextflow src/nf_modules/cutadapt/trimming_paired.nf \
  -c src/nf_modules/cutadapt/trimming_paired.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq" \
  -resume

./nextflow src/nf_modules/cutadapt/trimming_single.nf \
  -c src/nf_modules/cutadapt/trimming_single.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny*_S.fastq" \
  -resume

if [ -x "$(command -v singularity)" ]; then
./nextflow src/nf_modules/cutadapt/adaptor_removal_paired.nf \
  -c src/nf_modules/cutadapt/adaptor_removal_paired.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq" \
  -resume

./nextflow src/nf_modules/cutadapt/adaptor_removal_single.nf \
  -c src/nf_modules/cutadapt/adaptor_removal_single.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny*_S.fastq" \
  -resume

./nextflow src/nf_modules/cutadapt/trimming_paired.nf \
  -c src/nf_modules/cutadapt/trimming_paired.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq" \
  -resume

./nextflow src/nf_modules/cutadapt/trimming_single.nf \
  -c src/nf_modules/cutadapt/trimming_single.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny*_S.fastq" \
  -resume
fi
