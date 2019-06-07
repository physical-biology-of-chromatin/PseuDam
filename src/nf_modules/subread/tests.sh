./nextflow src/nf_modules/subread/subread.nf \
  -c src/nf_modules/subread/subread.config \
  -profile docker \
  --gtf "data/tiny_dataset/annot/tiny.gff" \
  --bam "data/tiny_dataset/map/tiny_v2.bam" \
  -resume

if [ -x "$(command -v singularity)" ]; then
./nextflow src/nf_modules/subread/subread.nf \
  -c src/nf_modules/subread/subread.config \
  -profile singularity \
  --gtf "data/tiny_dataset/annot/tiny.gff" \
  --bam "data/tiny_dataset/map/tiny_v2.bam" \
  -resume
fi
