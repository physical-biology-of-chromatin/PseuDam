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

