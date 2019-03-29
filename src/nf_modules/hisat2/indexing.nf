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
hisat2-build -p {task.cpus} ${fasta} ${fasta.baseName}.index
"""
}
