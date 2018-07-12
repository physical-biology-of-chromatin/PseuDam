params.gtf = "/media/manu/ManuDisque/gencode/gencode.v28.annotation.gtf"
params.genome = "/media/manu/ManuDisque/gencode/GRCh38.p12.genome.fa"

log.info "gtf file : ${params.gtf}"
log.info "genome fasta file : ${params.genome}"

Channel
  .fromPath( params.gtf )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.gtf}" }
  .set { gtf_file }
Channel
  .fromPath( params.genome )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.genome}" }
  .set { genome_file }

process create_annot {
  publishDir "results/ribowave/", mode: 'copy'

  input:
    file gtf from gtf_file
    file genome from genome_file

  output:
    file "*" into annot_file

  script:
"""
mkdir annotation
/Ribowave/scripts/create_annotation.sh -G ${gtf} -f ${genome}  -o annotation  -s /Ribowave/scripts
"""
}
