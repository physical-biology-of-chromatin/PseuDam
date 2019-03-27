cp data/tiny_dataset/map/tiny_v2.sort.bam data/tiny_dataset/map/tiny_v2_control.sort.bam
./nextflow src/nf_modules/music/peak_calling_single.nf \
  -c src/nf_modules/music/peak_calling_single.config \
  -profile docker \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" \
  --bam "data/tiny_dataset/map/*.sort.bam" \
  --index "data/tiny_dataset/map/*.sort.bam.bai*" \
  --read_size 50 --frag_size 300 \
  -resume

if [ -x "$(command -v singularity)" ]; then
./nextflow src/nf_modules/music/peak_calling_single.nf \
  -c src/nf_modules/music/peak_calling_single.config \
  -profile singularity \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" \
  --bam "data/tiny_dataset/map/*.sort.bam" \
  --index "data/tiny_dataset/map/*.sort.bam.bai*" \
  --read_size 50 --frag_size 300 \
  -resume
fi
