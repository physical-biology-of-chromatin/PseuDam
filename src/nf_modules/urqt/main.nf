version = "d62c1f8"
container_url = "lbmc/urqt:${version}"

trim_quality = "20"

process trimming_pairedend {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "${reads}"

  input:
  tuple val(pair_id), path(reads)

  output:
  tuple val(pair_id), path("*_trim_R{1,2}.fastq.gz"), emit: fastq
  path "*_report.txt", emit: report

  script:
"""
UrQt --t 20 --m ${task.cpus} --gz \
  --in ${reads[0]} --inpair ${reads[1]} \
  --out ${pair_id}_trim_R1.fastq.gz --outpair ${pair_id}_trim_R2.fastq.gz \
  > ${pair_id}_trimming_report.txt
"""
}

process trimming_singleend {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
  tuple val(file_id), path(reads)

  output:
  tuple val(file_id), path("*_trim.fastq.gz"), emit: fastq
  path "*_report.txt", emit: report

  script:
"""
UrQt --t 20 --m ${task.cpus} --gz \
  --in ${reads} \
  --out ${file_id}_trim.fastq.gz \
  > ${file_id}_trimming_report.txt
"""
}

