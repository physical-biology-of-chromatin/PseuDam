/*
* for single-end data
*/

params.fastq = "$baseDir/data/fastq/*.fastq"
params.index = "$baseDir/data/index/*.index*"

log.info "fastq files : ${params.fastq}"
log.info "index files : ${params.index}"

Channel
  .fromPath( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }
Channel
  .fromPath( params.index )
  .ifEmpty { error "Cannot find any index files matching: ${params.index}" }
  .set { index_files }

process mapping_fastq {
  tag "$reads.baseName"
  cpus 4
  publishDir "results/mapping/", mode: 'copy'

  input:
  file reads from fastq_files
  file index from index_files.toList()

  output:
  file "*" into count_files

  script:
"""
hisat2 -x ${file(file(index[0]).baseName).baseName} -U ${reads} -S ${reads.baseName}.sam -p ${task.cpus}
"""
}
