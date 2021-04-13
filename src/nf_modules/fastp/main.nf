version = "0.20.1"
container_url = "lbmc/fastp:${version}"

params.fastp_protocol = ""
params.fastp = ""
params.fastp_pairedend = ""
params.fastp_singleend = ""

workflow fastp {
  take:
    fastq

  main:
  switch(params.fastp_protocol) {
    case "accel_1splus":
      fastp_accel_1splus(fastq)
      fastp_accel_1splus.out.fastq.set{res_fastq}
      fastp_accel_1splus.out.report.set{res_report}
    break;
    default:
      fastp_default(fastq)
      fastp_default.out.fastq.set{res_fastq}
      fastp_default.out.report.set{res_report}
    break;
  }
  emit:
    fastq = res_fastq
    report = res_report
}

process fastp_default {
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
if (reads instanceof List)
"""
fastp --thread ${task.cpus} \
  --qualified_quality_phred 20 \
  --disable_length_filtering \
  --detect_adapter_for_pe \
  ${params.fastp} \
  --in1 ${reads[0]} \
  --in2 ${reads[1]} \
  --out1 ${pair_id}_R1_trim.fastq.gz \
  --out2 ${pair_id}_R2_trim.fastq.gz \
  --html ${pair_id}.html \
  --json ${pair_id}_fastp.json \
  --report_title ${pair_id}
"""
else
"""
fastp --thread ${task.cpus} \
  --qualified_quality_phred 20 \
  --disable_length_filtering \
  --detect_adapter_for_pe \
  ${params.fastp} \
  --in1 ${reads} \
  --out1 ${pair_id}_trim.fastq.gz \
  --html ${pair_id}.html \
  --json ${pair_id}_fastp.json \
  --report_title ${pair_id}
"""
}

process fastp_accel_1splus {
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
if (reads instanceof List)
"""
fastp --thread ${task.cpus} \
  --disable_quality_filtering \
  --disable_length_filtering \
  --disable_trim_poly_g \
  --stdout \
  --in1 ${reads[0]} \
  --in2 ${reads[1]} \
  --out1 ${pair_id}_R1_trim.fastq.gz \
  --out2 ${pair_id}_R2_trim.fastq.gz | \
  fastp --thread ${task.cpus} \
    --stdin \
    --interleaved_in \
    --trim_front1=10 \
    --trim_front2=10 \
    --qualified_quality_phred 20 \
    --disable_length_filtering \
    --detect_adapter_for_pe \
    ${params.fastp} \
    --html ${pair_id}.html \
    --json ${pair_id}_fastp.json \
    --report_title ${pair_id}
"""
else
"""
fastp --thread ${task.cpus} \
  --disable_quality_filtering \
  --disable_length_filtering \
  --disable_trim_poly_g \
  --stdout \
  --in1 ${reads[0]} \
  --out1 ${pair_id}_R1_trim.fastq.gz \
  fastp --thread ${task.cpus} \
    --stdin \
    --trim_front1=10 \
    --qualified_quality_phred 20 \
    --disable_length_filtering \
    --detect_adapter_for_pe \
    ${params.fastp} \
    --html ${pair_id}.html \
    --json ${pair_id}_fastp.json \
    --report_title ${pair_id}
"""
}

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
${params.fastp_pairedend} \
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
${params.fastp_singleend} \
--in1 ${reads} \
--out1 ${pair_id}_trim.fastq.gz \
--html ${pair_id}.html \
--json ${pair_id}_fastp.json \
--report_title ${pair_id}
"""
}
