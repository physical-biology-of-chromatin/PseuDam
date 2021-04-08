version = "2.18.11"
container_url = "lbmc/picard:${version}"

params.mark_duplicate = "VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true"
process mark_duplicate {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam)
  output:
    tuple val(file_id) , path("*.bam"), emit: bam
    path "*_report.txt", emit: report


  script:
"""
PicardCommandLine MarkDuplicates \
  ${params.mark_duplicate} \
  INPUT=${bam} \
  OUTPUT=${bam.baseName}_dedup.bam \
  METRICS_FILE=${bam.baseName}_picard_dedup_report.txt &> \
  picard_${bam.baseName}.log
"""
}

params.index_fasta = ""
process index_fasta {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(fasta)
  output:
    tuple val(file_id), path("*.dict"), emit: index

  script:
"""
PicardCommandLine CreateSequenceDictionary \
  ${params.index_fasta} \
  REFERENCE=${fasta} \
  OUTPUT=${fasta.baseName}.dict
"""
}

params.index_bam = ""
process index_bam {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam)
  output:
    tuple val(file_id), path("*"), emit: index

  script:
"""
PicardCommandLine BuildBamIndex \
  ${params.index_bam} \
  INPUT=${bam}
"""
}
