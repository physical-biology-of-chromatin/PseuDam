./nextflow src/nf_modules/multiqc/multiqc_paired.nf \
  -c src/nf_modules/multiqc/multiqc_paired.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq" \
  -resume

./nextflow src/nf_modules/multiqc/multiqc_single.nf \
  -c src/nf_modules/multiqc/multiqc_single.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_S.fastq" \
  -resume

if [ -x "$(command -v singularity)" ]; then
./nextflow src/nf_modules/multiqc/multiqc_paired.nf \
  -c src/nf_modules/multiqc/multiqc_paired.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq" \
  -resume

./nextflow src/nf_modules/multiqc/multiqc_single.nf \
  -c src/nf_modules/multiqc/multiqc_single.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_S.fastq" \
  -resume
fi
