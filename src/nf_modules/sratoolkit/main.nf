version = "2.8.2"
container_url = "lbmc/sratoolkit:${version}"

process index_fasta {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$sra"

  input:
    val sra

  output:
    tuple val(sra), path("*.fastq"), emit: fastq

  script:
"""
fastq-dump --split-files --gzip ${sra}
if [ -f ${sra}_1.fastq ]
then
  mv ${sra}_1.fastq ${sra}_R1.fastq
fi
if [ -f ${sra}_2.fastq ]
then
  mv ${sra}_2.fastq ${sra}_R2.fastq
fi
"""
}
