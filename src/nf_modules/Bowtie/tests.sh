./nextflow src/nf_modules/Bowtie/indexing.nf \
  -c src/nf_modules/Bowtie/indexing.config \
  -profile docker \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta"

./nextflow src/nf_modules/Bowtie/mapping_single.nf \
  -c src/nf_modules/Bowtie/mapping_single.config \
  -profile docker \
  --index "results/mapping/index/*.ebwt" \
  --fastq "data/tiny_dataset/fastq/tiny*_S.fastq"

./nextflow src/nf_modules/Bowtie/mapping_paired.nf \
  -c src/nf_modules/Bowtie/mapping_paired.config \
  -profile docker \
  --index "results/mapping/index/*.ebwt" \
  --fastq "data/tiny_dataset/fastq/tiny*_R{1,2}.fastq"

if [ -x "$(command -v singularity)" ]; then
./nextflow src/nf_modules/Bowtie/indexing.nf \
  -c src/nf_modules/Bowtie/indexing.config \
  -profile singularity \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta"

./nextflow src/nf_modules/Bowtie/mapping_single.nf \
  -c src/nf_modules/Bowtie/mapping_single.config \
  -profile singularity \
  --index "results/mapping/index/*.ebwt" \
  --fastq "data/tiny_dataset/fastq/tiny*_S.fastq"

./nextflow src/nf_modules/Bowtie/mapping_paired.nf \
  -c src/nf_modules/Bowtie/mapping_paired.config \
  -profile singularity \
  --index "results/mapping/index/*.ebwt" \
  --fastq "data/tiny_dataset/fastq/tiny*_R{1,2}.fastq"
fi
