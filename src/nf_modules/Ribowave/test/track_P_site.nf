params.bam = ""
params.exon = ""
params.genome = ""
params.jobname = ""
params.p_site = ""

log.info "bam file(s) : ${params.bam}"
log.info "exon file : ${params.exon}"
log.info "genome file : ${params.genome}"
log.info "job name : ${params.jobname}"
log.info "job name : ${params.p_site}"

Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.bam}" }
  .set { bam_files }
Channel
  .fromPath( params.exon )
  .ifEmpty { error "Cannot find any index files matching: ${params.exon}" }
  .set { exon_file }
Channel
  .fromPath( params.genome )
  .ifEmpty { error "Cannot find any index files matching: ${params.genome}" }
  .set { genome_file }
Channel
  .fromPath( params.p_site )
  .ifEmpty { error "Cannot find any index files matching: ${params.p_site}" }
  .set { p_site_file }

process determination_P_site {
  publishDir "results/ribowave/track_P_site", mode: 'copy'

  input:
  file bam from bam_files
  file exon from exon_file
  file genome from genome_file
  file p_site from p_site_file

  output:
  file "*" into det_p_site_channel

  script:
"""
/Ribowave/scripts/create_track_Ribo.sh -i ${bam} -G ${exon} -g ${genome} -P ${p_site} -o ./ -n ${params.jobname} -s /Ribowave/scripts
"""
}
