params.fastq = "$baseDir/data/fastq/*_{1,2}.fastq"
params.index = "$baseDir/data/index/*.index.*"

log.info "fastq files : ${params.fastq}"
log.info "index files : ${params.index}"

Channel
  .fromFilePairs( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }
Channel
  .fromPath( params.index )
  .ifEmpty { error "Cannot find any index files matching: ${params.index}" }
  .set { index_files }

process mapping_fastq {
  tag "$reads"
  //tag "$index.baseName"
  cpus 4
  publishDir "results/mapping/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files
  file index from index_files.toList()

  output:
  file "*" into counts_files

  script:
"""
hisat2 -x ${file(file(index[0]).baseName).baseName} -1 ${reads[0]} -2 ${reads[1]} -S ${pair_id}.sam -p ${task.cpus}
"""
}
