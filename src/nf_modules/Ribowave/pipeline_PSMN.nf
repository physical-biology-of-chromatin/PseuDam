/*
* Ribowave :
* Inputs : gtf genome files
* Inputs : bam file
* Inputs : genome size file
*/

/* 		PARAMETERS		 */

params.gtf = ""
params.genome = ""
params.bam = ""
params.genomesize = ""

log.info "gtf file : ${params.gtf}"
log.info "genome fasta file : ${params.genome}"
log.info "bam file(s) : ${params.bam}"
log.info "genomesize file : ${params.genomesize}"

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
bam_files.into {bam_deter_P_site ; bam_track_P_site ; bam_ribowave}
Channel
  .fromPath( params.genomesize )
  .ifEmpty { error "Cannot find any index files matching: ${params.genomesize}" }
  .set { genomesize_file }


/*		CREATE ANNOTATION		*/

process create_annot {
  publishDir "results/ribowave/annotation", mode: 'copy'

  input:
    file gtf from gtf_file
    file genome from genome_file

  output:
    file "*" into annot_file
    file "start_codon.bed" into start_codon_channel
    file "exons.gtf" into exon_gtf_channel
    file "final.ORFs" into finalORF_channel

  script:
"""
/Xnfs/lbmcdb/Ricci_team/shared_data/softwares/Ribowave/scripts/create_annotation.sh -G ${gtf} -f ${genome}  -o ./  -s /Xnfs/lbmcdb/Ricci_team/shared_data/softwares/Ribowave/scripts
"""
}


/*		P-site determination		*/

process determination_P_site {
  tag "$bam.baseName"
  publishDir "results/ribowave/", mode: 'copy'

  input:
  file bam from bam_deter_P_site
  file start from start_codon_channel

  output:
  file "*" into p_site_channel
  file "P-site/*1nt.txt" into psite1nt_channel

  script:
"""
/Xnfs/lbmcdb/Ricci_team/shared_data/softwares/Ribowave/scripts/P-site_determination.sh -i ${bam} -S ${start} -o ./ -n ${bam.baseName} -s /Xnfs/lbmcdb/Ricci_team/shared_data/softwares/Ribowave/scripts
"""
}

/*		P-site track		*/
process track_P_site {
  tag "$bam.baseName"
  publishDir "results/ribowave", mode: 'copy'

  input:
  file bam from bam_track_P_site
  file exon from exon_gtf_channel
  file genomesize from genomesize_file
  file p_site from psite1nt_channel

  output:
  file "*" into track_p_site_channel
  file "bedgraph/${bam.baseName}/final.psite" into psite_channel 
  

  script:
"""
/Xnfs/lbmcdb/Ricci_team/shared_data/softwares/Ribowave/scripts/create_track_Ribo.sh -i ${bam} -G ${exon} -g ${genomesize} -P ${p_site} -o ./ -n ${bam.baseName} -s /Xnfs/lbmcdb/Ricci_team/shared_data/softwares/Ribowave/scripts
"""
}
/*		ribowave Identifying translated ORF		*/

process ribowave_transORF {
  tag "$bam.baseName"
  publishDir "results/ribowave", mode: 'copy'

  input:
  file psite from psite_channel
  file bam from bam_ribowave
  file finalORF from finalORF_channel

  output:
  file "*" into ribowave_channel

  script:
"""
/Xnfs/lbmcdb/Ricci_team/shared_data/softwares/Ribowave/scripts/Ribowave -PD -a ${psite} -b ${finalORF} -o ./ -n ${bam.baseName} -s /Xnfs/lbmcdb/Ricci_team/shared_data/softwares/Ribowave/scripts -p 16
"""
}

