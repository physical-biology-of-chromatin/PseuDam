nextflow src/nf_modules/UrQt/tests/trimming_paired.nf \
  -c src/nf_modules/UrQt/urqt.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq"

nextflow src/nf_modules/UrQt/tests/trimming_single.nf \
  -c src/nf_modules/UrQt/urqt.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq"