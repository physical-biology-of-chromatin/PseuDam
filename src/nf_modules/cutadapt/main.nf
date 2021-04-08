version = "2.1"
container_url = "lbmc/cutadapt:${version}"

adapter_3_prim = "AGATCGGAAGAG"
adapter_5_prim = "CTCTTCCGATCT"
trim_quality = "20"

params.adaptor_removal = "-a ${adapter_3_prim} -g ${adapter_5_prim} -A ${adapter_3_prim} -G ${adapter_5_prim}"
process adaptor_removal {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$pair_id"

  input:
  tuple val(pair_id), path(reads)

  output:
  tuple val(pair_id), path("*_cut_R{1,2}.fastq.gz"), emit: fastq
  path "*_report.txt", emit: report

  script:
if (reads instanceof List)
  """
  cutadapt ${params.adaptor_removal} \
  -o ${pair_id}_cut_R1.fastq.gz -p ${pair_id}_cut_R2.fastq.gz \
  ${reads[0]} ${reads[1]} > ${pair_id}_report.txt
  """
else:
  """
  cutadapt ${params.adaptor_removal} \
  -o ${file_id}_cut.fastq.gz \
  ${reads} > ${file_id}_report.txt
  """
}

params.adaptor_removal_pairedend = "-a ${adapter_3_prim} -g ${adapter_5_prim} -A ${adapter_3_prim} -G ${adapter_5_prim}"
process adaptor_removal_pairedend {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$pair_id"

  input:
  tuple val(pair_id), path(reads)

  output:
  tuple val(pair_id), path("*_cut_R{1,2}.fastq.gz"), emit: fastq
  path "*_report.txt", emit: report

  script:
  """
  cutadapt ${params.adaptor_removal_pairedend} \
  -o ${pair_id}_cut_R1.fastq.gz -p ${pair_id}_cut_R2.fastq.gz \
  ${reads[0]} ${reads[1]} > ${pair_id}_report.txt
  """
}

params.adaptor_removal_singleend = "-a ${adapter_3_prim} -g ${adapter_5_prim}"
process adaptor_removal_singleend {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
  tuple val(file_id), path(reads)

  output:
  tuple val(file_id), path("*_cut.fastq.gz"), emit: fastq
  path "*_report.txt", emit: report

  script:
  """
  cutadapt ${params.adaptor_removal_singleend} \
  -o ${file_id}_cut.fastq.gz \
  ${reads} > ${file_id}_report.txt
  """
}

process trimming_pairedend {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$pair_id"

  input:
  tuple val(pair_id), path(reads)

  output:
  tuple val(pair_id), path("*_trim_R{1,2}.fastq.gz"), emit:fastq
  path "*_report.txt", emit: report

  script:
if (reads instanceof List)
  """
  cutadapt -q ${trim_quality},${trim_quality} \
  -o ${pair_id}_trim_R1.fastq.gz -p ${pair_id}_trim_R2.fastq.gz \
  ${reads[0]} ${reads[1]} > ${pair_id}_report.txt
  """
else
  """
  cutadapt -q ${trim_quality},${trim_quality} \
  -o ${file_id}_trim.fastq.gz \
  ${reads} > ${file_id}_report.txt
  """
}

process trimming_pairedend {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$pair_id"

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

