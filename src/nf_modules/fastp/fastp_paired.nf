params.fastq = "$baseDir/data/fastq/*_{1,2}.fastq"

log.info "fastq files : ${params.fastq}"

Channel
  .fromFilePairs( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

process fastp_fastq {
  tag "$pair_id"
  publishDir "results/fastq/fastp/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files

  output:
    file "*.{zip,html}" into fastp_report
    set pair_id, file "*.fastq.gz" fastq_trim_files

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
--report_title ${pair_id}
"""
}
