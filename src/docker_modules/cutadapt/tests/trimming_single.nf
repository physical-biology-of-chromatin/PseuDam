log.info "fastq files : ${params.fastq}"

Channel
  .fromPath( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

process trimming {
  tag "$reads.baseName"

  input:
  file reads from fastq_files

  output:
  file "*_trim.fastq.gz" into fastq_files_cut

  script:
  """
  cutadapt -q 20,20 \
  -o ${reads.baseName}_trim.fastq.gz \
  ${reads} > ${reads.baseName}_report.txt
  """
}

