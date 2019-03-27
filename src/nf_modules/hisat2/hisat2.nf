/*
* Hisat2 :
* Imputs : fastq files
* Imputs : fasta files
* Output : bam files
*/

/*                      fasta indexing                                     */
params.fasta = "$baseDir/data/bam/*.fasta"

log.info "fasta files : ${params.fasta}"

Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
  .set { fasta_file }

process index_fasta {
  tag "$fasta.baseName"
  publishDir "results/mapping/index/", mode: 'copy'

  input:
    file fasta from fasta_file

  output:
    file "*.index*" into index_files

  script:
"""
hisat2-build ${fasta} ${fasta.baseName}.index
"""
}

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

/*
* for paired-end data
*/

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

/*
* converting sam into bam
*/

/*                      sam to bam                                    */
params.sam = "$baseDir/data/bam/*.sam"

log.info "sam files : ${params.sam}"

Channel
  .fromPath( params.sam )
  .ifEmpty { error "Cannot find any sam files matching: ${params.sam}" }
  .set { sam_files }

process bam_converter {
  tag "$sam"
  cpus 4
  publishDir "results/mapping/bam/", mode: 'copy'

  input:
    file sam from sam_files

  output:
    file "*.bam" into bam_files

  script:
"""
samtools view -@ ${task.cpus} -bS ${sam} > ${sam.baseName}.bam
"""
}
