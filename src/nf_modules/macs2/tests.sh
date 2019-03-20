cp data/tiny_dataset/map/tiny_v2.bam data/tiny_dataset/map/tiny_v2_control.bam
./nextflow src/nf_modules/macs2/peak_calling.nf \
  -c src/nf_modules/macs2/peak_calling.config \
  -profile docker \
  -resume \
  --bam "data/tiny_dataset/map/tiny_v2*.bam" \
  --genome_size 129984 \
  --control_tag "control"
