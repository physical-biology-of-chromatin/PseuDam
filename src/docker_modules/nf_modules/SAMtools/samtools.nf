/*
* SAMtools :
* Imputs : bam files
* Output : bam files
*/

/*                      bams sorting                                     */
params.bam = "$baseDir/data/bam/*.bam"

log.info "bams files : ${params.bam}"

Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any bam files matching: ${params.bam}" }
  .set { bam_files }

process sort_bam {
  tag "$bam.baseName"
  cpus 4

  input:
    file bam from bam_files

  output:
    file "*_sorted.bam" into sorted_bam_files

  script:
"""
samtools sort -@ ${task.cpus} -O BAM -o ${bam.baseName}_sorted.bam ${bam}
"""
}

/*                      bams indexing                                     */

params.bam = "$baseDir/data/bam/*.bam"

log.info "bams files : ${params.bam}"

Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any bam files matching: ${params.bam}" }
  .set { bam_files }

process index_bam {
  tag "$bam.baseName"
  input:
    file bam from bam_files
  output:
    file "*bam*" into indexed_bam_file
  script:
"""
samtools index ${bam}
"""
}


/*                      bams spliting                                     */
params.bam = "$baseDir/data/bam/*.bam"

log.info "bams files : ${params.bam}"

Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any bam files matching: ${params.bam}" }
  .set { bam_files }

process split_bam {
  tag "$bam.baseName"
  cpus 2

  input:
    file bam from bam_files

  output:
    file "*_forward.bam*" into forward_bam_files
    file "*_reverse.bam*" into reverse_bam_files
  script:
"""
samtools view -hb -F 0x10 ${bam} > ${bam}_forward.bam &
samtools view -hb -f 0x10 ${bam} > ${bam}_reverse.bam
"""
}


/*                      bams filtering                                     */
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


