version = "2.25.0"
container_url = "lbmc/bedtools:${version}"

process fasta_from_bed {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${bed.baseName}"
  publishDir "results/fasta/", mode: 'copy'

  input:
  path fasta
  path bed

  output:
  path "*_extracted.fasta", emit: fasta

  script:
"""
bedtools getfasta -name \
-fi ${fasta} -bed ${bed} -fo ${bed.baseName}_extracted.fasta
"""
}
