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

