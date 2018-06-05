params.bam = "$baseDir/data/bam/*.bam"
params.bed = "$baseDir/data/bam/*.bed"

log.info "bams files : ${params.bam}"
log.info "bed file : ${params.bed}"

Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any bam files matching: ${params.bam}" }
  .set { bam_files }
Channel
  .fromPath( params.bed )
  .ifEmpty { error "Cannot find any bed file matching: ${params.bed}" }
  .set { bed_files }

process filter_bam {
  tag "$bam.baseName"
  cpus 4

  input:
    file bam from bam_files
    file bed from bed_files

  output:
    file "*_filtered.bam*" into filtered_bam_files
  script:
"""
samtools view -@ ${task.cpus} -hb ${bam} -L ${bed} > ${bam.baseName}_filtered.bam
"""
}


