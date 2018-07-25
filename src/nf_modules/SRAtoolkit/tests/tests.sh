nextflow src/nf_modules/SRAtoolkit/tests/fastqdump.nf \
  -c src/nf_modules/SRAtoolkit/sratoolkit.config \
  -profile docker \
  --list_srr "src/nf_modules/SRAtoolkit/tests/list-srr.txt"
