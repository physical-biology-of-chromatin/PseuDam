./nextflow src/nf_modules/bedtools/fasta_from_bed.nf \
  -c src/nf_modules/bedtools/fasta_from_bed.config \
  -profile docker \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" \
  --bed "data/tiny_dataset/annot/tiny.bed" \
  -resume

if [ -x "$(command -v singularity)" ]; then
./nextflow src/nf_modules/bedtools/fasta_from_bed.nf \
  -c src/nf_modules/bedtools/fasta_from_bed.config \
  -profile singularity \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" \
  --bed "data/tiny_dataset/annot/tiny.bed" \
  -resume
fi