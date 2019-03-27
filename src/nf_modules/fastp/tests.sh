./nextflow src/nf_modules/fastp/fastp_paired.nf \
  -c src/nf_modules/fastp/fastp_paired.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq" \
  -resume

./nextflow src/nf_modules/fastp/fastp_single.nf \
  -c src/nf_modules/fastp/fastp_single.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_S.fastq" \
  -resume

if [ -x "$(command -v singularity)" ]; then
./nextflow src/nf_modules/fastp/fastp_paired.nf \
  -c src/nf_modules/fastp/fastp_paired.config \
  -profile singularity \
  --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq" \
  -resume

./nextflow src/nf_modules/fastp/fastp_single.nf \
  -c src/nf_modules/fastp/fastp_single.config \
  -profile singularity \
  --fastq "data/tiny_dataset/fastq/tiny_S.fastq" \
  -resume
fi
