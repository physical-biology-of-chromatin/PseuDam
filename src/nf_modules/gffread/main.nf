version = "0.12.2"
container_url = "lbmc/gffread:${version}"

params.gffread = ""
params.gffread_out = ""
process gffread {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_prefix"
  if (params.gffread_out != "") {
    publishDir "results/${params.gffread_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(gtf)
  tuple val(fasta_id), path(fasta)

  output:
    tuple val(fasta_id), path("${file_prefix}.fasta"), emit: fasta

  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
  """
  gffread ${gtf} -g ${fasta} -o ${file_prefix}.fasta
  """
}