./nextflow src/nf_modules/Bowtie/tests/index.nf \
  -c src/nf_modules/Bowtie/bowtie.config \
  -profile docker \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta"

./nextflow src/nf_modules/Bowtie/tests/mapping_single.nf \
  -c src/nf_modules/Bowtie/bowtie.config \
  -profile docker \
  --index "results/mapping/index/*.ebwt" \
  --fastq "data/tiny_dataset/fastq/tiny*_S.fastq"

./nextflow src/nf_modules/Bowtie/tests/mapping_paired.nf \
  -c src/nf_modules/Bowtie/bowtie.config \
  -profile docker \
  --index "results/mapping/index/*.ebwt" \
  --fastq "data/tiny_dataset/fastq/tiny*_R{1,2}.fastq"

