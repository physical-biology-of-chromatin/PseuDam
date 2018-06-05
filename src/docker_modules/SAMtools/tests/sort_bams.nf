params.bam = "$baseDir/data/bam/*.bam"

log.info "bams files : ${params.bam}"

Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any bam files matching: ${params.bam}" }
  .set { bam_files }

process sort_bam {
  tag "$bams.baseName"
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

