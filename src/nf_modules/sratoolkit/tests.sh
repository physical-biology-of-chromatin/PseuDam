./nextflow src/nf_modules/sratoolkit/fastqdump.nf \
  -c src/nf_modules/sratoolkit/fastqdump.config \
  -profile docker \
  --list_srr "src/nf_modules/sratoolkit/list-srr.txt" \
  -resume

if [ -x "$(command -v singularity)" ]; then
./nextflow src/nf_modules/sratoolkit/fastqdump.nf \
  -c src/nf_modules/sratoolkit/fastqdump.config \
  -profile docker \
  --list_srr "src/nf_modules/sratoolkit/list-srr.txt" \
  -resume
fi
