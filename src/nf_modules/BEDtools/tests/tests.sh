nextflow src/nf_modules/BEDtools/tests/fasta_from_bed.nf \
  -c src/nf_modules/BEDtools/bedtools.config \
  -profile docker \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" \
  --bed "data/tiny_dataset/annot/tiny.bed" \
