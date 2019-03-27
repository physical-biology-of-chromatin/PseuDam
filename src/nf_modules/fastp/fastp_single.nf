params.fastq = "$baseDir/data/fastq/*.fastq"

log.info "fastq files : ${params.fastq}"

Channel
  .fromPath( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { fastq_files }

process fastp_fastq {
  tag "$file_id"
  publishDir "results/fastq/fastp/", mode: 'copy'

  input:
  set file_id, file(reads) from fastq_files

  output:
    file "*.{zip,html}" into fastp_report
    set file_id, file "*.fastq.gz" fastq_trim_files

  script:
"""
fastp --thread ${task.cpus} \
--qualified_quality_phred 20 \
--disable_length_filtering \
--in1 ${reads} \
--out1 ${file_id}_R1_trim.fastq.gz \
--html ${file_id}.html \
--report_title ${file_id}
"""
}
