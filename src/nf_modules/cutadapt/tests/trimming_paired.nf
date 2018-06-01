log.info "fastq files : ${params.fastq}"

Channel
  .fromFilePairs( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

process trimming {
  tag "$pair_id"

  input:
  set pair_id, file(reads) from fastq_files

  output:
  file "*_trim_R{1,2}.fastq.gz" into fastq_files_cut

  script:
  """
  cutadapt -q 20,20 \
  -o ${pair_id}_trim_R1.fastq.gz -p ${pair_id}_trim_R2.fastq.gz \
  ${reads[0]} ${reads[1]} > ${pair_id}_report.txt
  """
}

