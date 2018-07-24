params.bam = ""
params.start = ""

log.info "bam file(s) : ${params.bam}"
log.info "start_codon file : ${params.start}"

Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.bam}" }
  .set { bam_files }
Channel
  .fromPath( params.start )
  .ifEmpty { error "Cannot find any index files matching: ${params.start}" }
  .set { start_file }

process determination_P_site {
  tag "$bam.baseName"
  publishDir "results/ribowave", mode: 'copy'

  input:
  file bam from bam_files
  file start from start_file

  output:
  file "*" into p_site_channel
  file "*psite1nt.txt" into psite1nt_channel

  script:
"""
/Ribowave/scripts/P-site_determination.sh -i ${bam} -S ${start} -o ./ -n ${bam.baseName} -s /Ribowave/scripts
"""
}
