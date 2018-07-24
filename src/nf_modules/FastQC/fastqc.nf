/*
* fastqc :
* Imputs : fastq files
* Output : pdf files
*/

/*                      fastQC                                     */
params.fastq = ""

log.info "fastq files : ${params.fastq}"

Channel
  .fromPath( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

process fastqc {
  tag "$fastq.baseName"
  publishDir "results/fasqc/${fastq.baseName}/", mode: 'copy'

  input:
    file fastq from fastq_files

  output:
    file "*" into fastqc_repport

  script:
"""
fastqc -o ./ --noextract -f fastq ${fastq}
"""
}


