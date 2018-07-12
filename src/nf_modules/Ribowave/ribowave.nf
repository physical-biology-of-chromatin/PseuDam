/*
* Ribowave :
* Inputs : fastq files
* Inputs : fasta files
* Output : bam files
*/

/*                      Create annotation                                     */
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
  publishDir "results/ribowave/annotation", mode: 'copy'

  input:
    file gtf from gtf_file
    file genome from genome_file

  output:
    file "*" into annot_file

  script:
"""
/Ribowave/scripts/create_annotation.sh -G ${gtf} -f ${genome}  -o ./  -s /Ribowave/scripts
"""
}


/*
* P-site determination
*/

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

process determination_P_site {
  publishDir "results/ribowave/", mode: 'copy'

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

/*
* P-site track
*/

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

process track_P_site {
  publishDir "results/ribowave", mode: 'copy'

  input:
  file bam from bam_files
  file exon from exon_file
  file genome from genome_file
  file p_site from p_site_file

  output:
  file "*" into track_p_site_channel

  script:
"""
/Ribowave/scripts/create_track_Ribo.sh -i ${bam} -G ${exon} -g ${genome} -P ${p_site} -o ./ -n ${params.jobname} -s /Ribowave/scripts
"""
}

/*
* ribowave Identifying translated ORF
*/

params.psite = ""
params.finalORF = ""
params.outputdir = ""
params.jobname = ""

log.info "psite file : ${params.psite}"
log.info "finalORF file : ${params.finalORF}"
log.info "job name : ${params.jobname}"
log.info "output dir : ${params.outputdir}"

Channel
  .fromPath( params.psite )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.psite}" }
  .set { psite_file }
Channel
  .fromPath( params.finalORF )
  .ifEmpty { error "Cannot find any index files matching: ${params.finalORF}" }
  .set { finalORF_file }

process ribowave_transORF {
  publishDir "results/ribowave", mode: 'copy'

  input:
  file psite from psite_file
  file finalORF from finalORF_file

  output:
  file "*" into ribowave_channel

  script:
"""
/Ribowave/scripts/Ribowave -PD -a ${psite} -b ${finalORF} -o ${params.outputdir} -n ${params.jobname} -s /Ribowave/scripts -p 8
"""
}

