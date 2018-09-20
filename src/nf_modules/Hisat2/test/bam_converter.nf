/*
* SAMtools :
* Imputs : sam files
* Output : bam files
*/

/*                      sam to bam                                    */
params.sam = "$baseDir/data/bam/*.sam"

log.info "sam files : ${params.sam}"

Channel
  .fromPath( params.sam )
  .ifEmpty { error "Cannot find any sam files matching: ${params.sam}" }
  .set { sam_files }

process bam_converter {
  tag "$sam"
  cpus 4
  publishDir "results/mapping/bam/", mode: 'copy'

  input:
    file sam from sam_files

  output:
    file "*.bam" into bam_files

  script:
"""
samtools view -@ ${task.cpus} -bS ${sam} > ${sam.baseName}.bam
"""
}
