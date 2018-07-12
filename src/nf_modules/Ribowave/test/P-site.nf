params.bam = ""
params.start = ""
params.jobname = ""

log.info "bam file(s) : ${params.bam}"
log.info "start_codon file : ${params.start}"
log.info "job name : ${params.jobname}"

Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.bam}" }
  .set { bam_files }
Channel
  .fromPath( params.start )
  .ifEmpty { error "Cannot find any index files matching: ${params.start}" }
  .set { start_file }

process p_site {
  publishDir "results/ribowave/P-site", mode: 'copy'

  input:
  file bam from bam_files
  file start from start_file

  output:
  file "*" into p_site_channel

  script:
"""
/Ribowave/scripts/P-site_determination.sh -i ${bam} -S ${start} -o ./ -n ${params.jobname} -s /Ribowave/scripts
"""
}
