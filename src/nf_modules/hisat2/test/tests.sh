nextflow src/nf_modules/Hisat2/test/index.nf \
  -c src/nf_modules/Hisat2/hisat2.config \
  -profile docker \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta"

nextflow src/nf_modules/Hisat2/test/mapping_paired.nf \
  -c src/nf_modules/Hisat2/hisat2.config \
  -profile docker \
  --index "results/mapping/index/tiny_v2.index*" \
  --fastq "data/tiny_dataset/fastq/tiny*_R{1,2}.fastq"

nextflow src/nf_modules/Hisat2/test/mapping_single.nf \
  -c src/nf_modules/Hisat2/hisat2.config \
  -profile docker \
  --index "results/mapping/index/tiny_v2.index*" \
  --fastq "data/tiny_dataset/fastq/tiny*_S.fastq"

nextflow src/nf_modules/Hisat2/test/bam_converter.nf \
  -c src/nf_modules/Hisat2/hisat2.config \
  -profile docker \
  --sam "results/mapping/*.sam" \
