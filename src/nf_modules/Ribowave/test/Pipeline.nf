/*
* Ribowave :
* Inputs : fastq files
* Inputs : fasta files
* Output : bam files
*/

params.gtf = "/media/manu/ManuDisque/gencode/gencode.v28.annotation.gtf"
params.genome = "/media/manu/ManuDisque/gencode/GRCh38.p12.genome.fa"
params.bam = ""
params.jobname = ""

log.info "gtf file : ${params.gtf}"
log.info "genome fasta file : ${params.genome}"
log.info "bam file(s) : ${params.bam}"
log.info "job name : ${params.jobname}"

Channel
  .fromPath( params.gtf )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.gtf}" }
  .set { gtf_file }
Channel
  .fromPath( params.genome )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.genome}" }
  .set { genome_file }
Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.bam}" }
  .set { bam_files }

process create_annot {
  publishDir "results/ribowave/gtf27/annot", mode: 'copy'

  input:
    file gtf from gtf_file
    file genome from genome_file

  output:
    file "*" into annot_file_save
    file "start_codon.bed" into annot_file

  script:
"""
/Ribowave/scripts/create_annotation.sh -G ${gtf} -f ${genome} -o ./ -s /Ribowave/scripts
"""
}

process determination_P_site {
  publishDir "results/ribowave/gtf27/deter_P_site", mode: 'copy'

  input:
  file bam from bam_files
  file start from annot_file

  output:
  file "*" into p_site_channel

  script:
"""
/Ribowave/scripts/P-site_determination.sh -i ${bam} -S ${start} -o ./ -n ${params.jobname} -s /Ribowave/scripts
"""
}

