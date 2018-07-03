nextflow src/nf_modules/sra-tools/tests/fastqdump.nf \
  -c src/nf_modules/sra-tools/sra-tools.config \
  -profile docker \
  --list_srr "src/nf_modules/sra-tools/tests/list-srr.txt"
