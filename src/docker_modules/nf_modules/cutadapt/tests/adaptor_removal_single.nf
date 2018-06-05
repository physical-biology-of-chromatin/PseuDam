log.info "fastq files : ${params.fastq}"

Channel
  .fromPath( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

process adaptor_removal {
  tag "$reads.baseName"

  input:
  file reads from fastq_files

  output:
  file "*_cut.fastq.gz" into fastq_files_cut

  script:
  """
  cutadapt -a AGATCGGAAGAG -g CTCTTCCGATCT\
  -o ${reads.baseName}_cut.fastq.gz \
  ${reads} > ${reads.baseName}_report.txt
  """
}

