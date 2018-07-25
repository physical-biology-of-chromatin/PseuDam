/*
* fastqc :
* Imputs : fastq files
* Output : pdf files
*/

/*                      fastQC                                     */


/*
* for single-end data
*/

params.fastq = "$baseDir/data/fastq/*.fastq"

log.info "fastq files : ${params.fastq}"

Channel
  .fromPath( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

process fastqc_fastq {
  tag "$fastq.baseName"
  publishDir "results/fastq/fastqc/", mode: 'copy'
  cpus = 1

  input:
    file fastq from fastq_files

  output:
    file "*.{zip,html}" into fastqc_repport

  script:
"""
fastqc --quiet --threads ${task.cpus} --outdir -f fastq ./ ${fastq}
"""
}


/*
* for paired-end data
*/

params.fastq = "$baseDir/data/fastq/*_{1,2}.fastq"

log.info "fastq files : ${params.fastq}"

Channel
  .fromFilePairs( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

process fastqc_fastq {
  tag "$pair_id"
  publishDir "results/fastq/fastqc/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files

  output:
    file "*.{zip,html}" into fastqc_repport

  script:
"""
fastqc --quiet --threads ${task.cpus} --outdir -f fastq ./ \
${fastq[0]} ${fastq[1]}
"""
}

