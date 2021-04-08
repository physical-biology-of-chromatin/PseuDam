version = "0.11.5"
container_url = "lbmc/fastqc:${version}"

params.fastqc_fastq = ""
process fastqc_fastq {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$pair_id"

  input:
  tuple val(pair_id), path(reads)

  output:
  path "*.{zip,html}", emit: report

  script:
if (reads instanceof List)
"""
fastqc --quiet --threads ${task.cpus} --format fastq --outdir ./ \
  ${params.fastqc_fastq} \
  ${reads[0]} ${reads[1]}
"""
else
"""
  fastqc --quiet --threads ${task.cpus} --format fastq --outdir ./ ${reads}
"""
}

params.fastqc_fastq_pairedend = ""
process fastqc_fastq_pairedend {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$pair_id"

  input:
  tuple val(pair_id), path(reads)

  output:
  path "*.{zip,html}", emit: report

  script:
"""
fastqc --quiet --threads ${task.cpus} --format fastq --outdir ./ \
  ${params.fastqc_fastq_pairedend} \
  ${reads[0]} ${reads[1]}
"""
}

params.fastqc_fastq_singleend = ""
process fastqc_fastq_singleend {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
  tuple val(file_id), path(reads)

  output:
    path "*.{zip,html}", emit: report

  script:
"""
  fastqc --quiet --threads ${task.cpus} ${params.fastqc_fastq_singleend} --format fastq --outdir ./ ${reads}
"""
}

