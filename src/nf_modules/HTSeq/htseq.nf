/*
* htseq :
* Imputs : sorted bams files
* Imputs : gtf
* Output : counts files
*/
/*                      quality trimming                                     */

params.bam = "$baseDir/data/bam/*.bam"
params.gtf = "$baseDir/data/annotation/*.gtf"

log.info "bam files : ${params.bam}"
log.info "gtf files : ${params.gtf}"

Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.bam}" }
  .set { bam_files }
Channel
  .fromPath( params.gtf )
  .ifEmpty { error "Cannot find any gtf file matching: ${params.gtf}" }
  .set { gtf_file }

process counting {
  tag "$bam.baseName"
  publishDir "results/quantification/", mode: 'copy'

  input:
  file bam from bam_files
  file gtf from gtf_file

  output:
  file "*.count" into count_files

  script:
"""
htseq-count -r pos --mode=intersection-nonempty -a 10 -s no -t exon -i gene_id \
--format=bam ${bam} ${gtf} > ${bam.baseName}.count
"""
}
