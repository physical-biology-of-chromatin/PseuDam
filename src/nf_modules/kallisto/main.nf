version = "0.44.0"
container_url = "lbmc/kallisto:${version}"

params.index_fasta = "-k 31 --make-unique"
params.index_fasta_out = ""
process index_fasta {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.index_fasta_out != "") {
    publishDir "results/${params.index_fasta_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta)

  output:
    tuple val(file_id), path("*.index*"), emit: index
    tuple val(file_id), path("*_report.txt"), emit: report

  script:
"""
kallisto index ${params.index_fasta} -i ${fasta.baseName}.index ${fasta} \
2> ${fasta.baseName}_kallisto_index_report.txt
"""
}

params.mapping_fastq = "--bias --bootstrap-samples 100"
params.mapping_fastq_out = ""
process mapping_fastq {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$pair_id"
  if (params.mapping_fastq_out != "") {
    publishDir "results/${params.mapping_fastq_out}", mode: 'copy'
  }

  input:
  tuple val(index_id), path(index)
  tuple val(file_id), path(reads)

  output:
  tuple val(file_id), path("${pair_id}"), emit: counts
  tuple val(file_id), path("*_report.txt"), emit: report

  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }

  if (reads.size() == 2)
  """
  mkdir ${file_prefix}
  kallisto quant -i ${index} -t ${task.cpus} \
  ${params.mapping_fastq} -o ${file_prefix} \
  ${reads[0]} ${reads[1]} &> ${file_prefix}_kallisto_mapping_report.txt
  """
  else
  """
  mkdir ${file_prefix}
  kallisto quant -i ${index} -t ${task.cpus} --single \
  ${params.mapping_fastq} -o ${file_prefix} \
  ${reads} &> ${file_prefix}_kallisto_mapping_report.txt
  """
}
