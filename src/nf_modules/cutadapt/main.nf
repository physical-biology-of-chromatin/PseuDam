version = "2.1"
container_url = "lbmc/cutadapt:${version}"

adapter_3_prim = "AGATCGGAAGAG"
adapter_5_prim = "CTCTTCCGATCT"
trim_quality = "20"


process adaptor_removal_paired {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$pair_id"
  publishDir "results/fastq/adaptor_removal/", mode: 'copy'

  input:
  tuple val(pair_id), path(reads)

  output:
  tuple val(pair_id), path("*_cut_R{1,2}.fastq.gz"), emit: fastq
  path "*_report.txt", emit: report

  script:
  """
  cutadapt -a ${adapter_3_prim} -g ${adapter_5_prim} -A ${adapter_3_prim} -G ${adapter_5_prim} \
  -o ${pair_id}_cut_R1.fastq.gz -p ${pair_id}_cut_R2.fastq.gz \
  ${reads[0]} ${reads[1]} > ${pair_id}_report.txt
  """
}

process adaptor_removal_singleend {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  publishDir "results/fastq/adaptor_removal/", mode: 'copy'

  input:
  tuple val(file_id), path(reads)

  output:
  tuple val(file_id), path("*_cut.fastq.gz"), emit: fastq
  path "*_report.txt", emit: report

  script:
  """
  cutadapt -a ${adapter_3_prim} -g ${adapter_5_prim} \
  -o ${file_id}_cut.fastq.gz \
  ${reads} > ${file_id}_report.txt
  """
}

process trimming_pairedend {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$pair_id"
  publishDir "results/fastq/trimming/", mode: 'copy'

  input:
  tuple val(pair_id), path(reads)

  output:
  tuple val(pair_id), path("*_trim_R{1,2}.fastq.gz"), emit:fastq
  path "*_report.txt", emit: report

  script:
  """
  cutadapt -q ${trim_quality},${trim_quality} \
  -o ${pair_id}_trim_R1.fastq.gz -p ${pair_id}_trim_R2.fastq.gz \
  ${reads[0]} ${reads[1]} > ${pair_id}_report.txt
  """
}

process trimming_singleend {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
  tuple val(file_id), path(reads)

  output:
  tuple val(file_id), path("*_trim.fastq.gz"), emit: fastq
  path "*_report.txt", emit: report

  script:
  """
  cutadapt -q ${trim_quality},${trim_quality} \
  -o ${file_id}_trim.fastq.gz \
  ${reads} > ${file_id}_report.txt
  """
}

