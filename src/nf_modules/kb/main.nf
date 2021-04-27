version = "0.26.0"
container_url = "lbmc/kb:${version}"

params.kb_protocol = "10x_v3"

params.index_fasta = ""
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
    tuple val(gtf_id), path(gtf)
    tuple val(t2g_id), path(transcript_to_gene)

  output:
    tuple val(file_id), path("*.idx"), emit: index
    tuple val(t2g_id), path("${transcript_to_gene}"), emit: t2g
    tuple val(file_id), path("*_report.txt"), emit: report

  script:
"""
kb ref \
  -i ${fasta.simpleName}.idx \
  -g ${transcript_to_gene} \
  ${params.index_fasta} \
  -f1 ${fasta} ${gtf} > ${fasta.simpleName}_kb_index_report.txt
"""
}

params.count = ""
params.count_out = ""
workflow count {
  take:
    index
    fastq
    transcript_to_gene
    whitelist

  main:
  whitelist
    .ifEmpty(["NO WHITELIST", 0])
    .set{ whitelist_optional }
  switch(params.kb_protocol) {
    case "marsseq":
      kb_marseq(index, fastq, transcript_to_gene, whitelist_optional)
      kb_marseq.out.counts.set{res_counts}
      kb_marseq.out.report.set{res_report}
    break;
    default:
      kb_default(index, fastq, transcript_to_gene, whitelist_optional)
      kb_default.out.counts.set{res_counts}
      kb_default.out.report.set{res_report}
    break;
  }
  emit:
    counts = res_counts
    report = res_report
}

process kb_default {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_prefix"
  if (params.kb_out != "") {
    publishDir "results/${params.kb_out}", mode: 'copy'
  }

  input:
  tuple val(index_id), path(index)
  tuple val(file_id), path(reads)
  tuple val(t2g_id), path(transcript_to_gene)
  tuple val(whitelist_id), path(whitelist)

  output:
  tuple val(file_id), path("${file_prefix}"), emit: counts
  tuple val(file_id), path("*_report.txt"), emit: report

  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
  def whitelist_param = ""
  if (whitelist_id != "NO WHITELIST"){
    whitelist_param = "-w ${white_list}"
  }

  if (reads.size() == 2)
  """
  mkdir ${file_prefix}
  kb count  -t ${task.cpus} \
    -m ${task.memory} \
    -i ${index} \
    -g ${transcript_to_gene} \
    ${whitelist_param} \
    -x 10XV3
    ${params.mapping_fastq} \
    -o ${file_prefix} \
    ${reads[0]} ${reads[1]} &> ${file_prefix}_kb_mapping_report.txt
  """
}

process kb_marseq {
  // With the MARS-Seq protocol, we have:
  // on the read 1: 4 nt of bc plate
  // on the read 2: 6 nt of bc cell, and 8 nt of UMI
  // this process expect that the bc plate is removed from the read 1
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_prefix"
  if (params.kb_out != "") {
    publishDir "results/${params.kb_out}", mode: 'copy'
  }

  input:
  tuple val(index_id), path(index)
  tuple val(file_id), path(reads)
  tuple val(t2g_id), path(transcript_to_gene)
  tuple val(whitelist_id), path(whitelist)

  output:
  tuple val(file_id), path("${file_prefix}"), emit: counts
  tuple val(file_id), path("*_report.txt"), emit: report

  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
  def whitelist_param = ""
  if (whitelist_id != "NO WHITELIST"){
    whitelist_param = "-w ${white_list}"
  }

  if (reads.size() == 2)
  """
  mkdir ${file_prefix}
  kb count  -t ${task.cpus} \
    -m ${task.memory} \
    -i ${index} \
    -g ${transcript_to_gene} \
    ${whitelist_param} \
    -x 1,0,6:1,6,14:1,14,0
    ${params.mapping_fastq} \
    -o ${file_prefix} \
    ${reads[0]} ${reads[1]} &> ${file_prefix}_kb_mapping_report.txt
  """
}