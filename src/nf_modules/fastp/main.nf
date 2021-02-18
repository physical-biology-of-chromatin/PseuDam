version = "0.20.1"
container_url = "lbmc/fastp:${version}"


process fastp_pairedend {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$pair_id"
  publishDir "results/QC/fastp/", mode: 'copy', pattern: "*.html"

  input:
  tuple val(pair_id), path(reads)

  output:
    tuple val(pair_id), path("*.fastq.gz"), emit: fastq
    tuple val(pair_id), path("*.html"), emit: html
    tuple val(pair_id), path("*.json"), emit: report

  script:
"""
fastp --thread ${task.cpus} \
--qualified_quality_phred 20 \
--disable_length_filtering \
--detect_adapter_for_pe \
--in1 ${reads[0]} \
--in2 ${reads[1]} \
--out1 ${pair_id}_R1_trim.fastq.gz \
--out2 ${pair_id}_R2_trim.fastq.gz \
--html ${pair_id}.html \
--json ${pair_id}_fastp.json \
--report_title ${pair_id}
"""
}

process fastp_singleend {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$pair_id"
  publishDir "results/QC/fastp/", mode: 'copy', pattern: "*.html"

  input:
  tuple val(pair_id), path(reads)

  output:
    tuple val(pair_id), path("*.fastq.gz"), emit: fastq
    tuple val(pair_id), path("*.html"), emit: html
    tuple val(pair_id), path("*.json"), emit: report

  script:
"""
fastp --thread ${task.cpus} \
--qualified_quality_phred 20 \
--disable_length_filtering \
--detect_adapter_for_pe \
--in1 ${reads} \
--out1 ${pair_id}_trim.fastq.gz \
--html ${pair_id}.html \
--json ${pair_id}_fastp.json \
--report_title ${pair_id}
"""
}
