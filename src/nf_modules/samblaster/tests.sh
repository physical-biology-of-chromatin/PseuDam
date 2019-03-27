./nextflow src/nf_modules/samblaster/dedup_sams.nf \
  -c src/nf_modules/samblaster/dedup_sams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.bam" \
  -resume

if [ -x "$(command -v singularity)" ]; then
./nextflow src/nf_modules/samblaster/dedup_sams.nf \
  -c src/nf_modules/samblaster/dedup_sams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.bam" \
  -resume
fi
