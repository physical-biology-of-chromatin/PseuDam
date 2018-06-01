params.fastq = "$baseDir/data/fastq/*.fastq"

log.info "fastq files : ${params.fastq}"

Channel
  .fromPath( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

process trimming {
  tag "$reads.baseName"
  cpus 4

  input:
  file reads from fastq_files

  output:
  file "*_trim.fastq.gz" into fastq_files_trim

  script:
  """
  UrQt --t 20 --m ${task.cpus} --gz \
  --in ${reads} \
  --out ${reads.baseName}_trim.fastq.gz \
  > ${reads.baseName}_trimming_report.txt
  """
}

