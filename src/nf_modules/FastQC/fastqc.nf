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
  publishDir "results/fastqc/", mode: 'copy'

  input:
    file fastq from fastq_files

  output:
    file "*.htlm" into fastqc_repport

  script:
"""
fastqc -o ./ --noextract -f fastq ${fastq}
"""
}


