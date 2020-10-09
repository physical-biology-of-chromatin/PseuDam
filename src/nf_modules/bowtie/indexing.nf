/*                      fasta indexing                                     */

params.fasta = "$baseDir/data/bam/*.fasta"

log.info "fasta files : ${params.fasta}"

Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any bam files matching: ${params.fasta}" }
  .set { fasta_file }

process index_fasta {
  tag "$fasta.baseName"
  publishDir "results/mapping/index/", mode: 'copy'
  label "bowtie"

  input:
    file fasta from fasta_file

  output:
    file "*.index*" into index_files
    file "*_report.txt" into indexing_report

  script:
"""
bowtie-build --threads ${task.cpus} -f ${fasta} ${fasta.baseName}.index &> ${fasta.baseName}_bowtie_report.txt

if grep -q "Error" ${fasta.baseName}_bowtie_report.txt; then
  exit 1
fi
"""
}

