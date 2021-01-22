version = "2.18.11"
container_url = "lbmc/picard:${version}"

process mark_duplicate {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  publishDir "results/mapping/ddup/", mode: 'copy'

  input:
    tuple val(file_id), path(bam)
  output:
    tuple val(file_id) , path("*.bam"), emit: bam
    path "*_report.txt", emit: report


  script:
"""
PicardCommandLine MarkDuplicates \
  VALIDATION_STRINGENCY=LENIENT \
  REMOVE_DUPLICATES=true \
  INPUT=${bams[0]} \
  OUTPUT=${bams[0].baseName}_dedup.bam \
  METRICS_FILE=${bams[0].baseName}_picard_dedup_report.txt &> \
  picard_${bams[0].baseName}.log
"""
}
