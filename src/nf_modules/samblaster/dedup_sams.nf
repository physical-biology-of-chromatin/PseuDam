params.sam = "$baseDir/data/sam/*.bam"

log.info "bam files : ${params.bam}"

Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any bam files matching: ${params.bam}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { bam_files }

process dedup_sam {
  tag "$file_id"
  cpus 4

  input:
    set file_id, file(bam) from bam_files

  output:
    set file_id, "*_dedup.bam*" into dedup_bam_files
  script:
"""
samtools view -h ${bam} | \
samblaster --addMateTags | \
samtools view -Sb - > ${file_id}_dedub.bam
"""
}


