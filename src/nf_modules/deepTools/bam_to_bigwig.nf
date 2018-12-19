params.bam = "$baseDir/data/bam/*.bam"

log.info "bams files : ${params.bam}"

Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any bam files matching: ${params.bam}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { bam_files }

bam_files.into{
  bam_files_index;
  bam_files_bigwig
  }

process index_bam {
  tag "$file_id"
  cpus 4

  input:
    set file_id, file(bam) from bam_files_index

  output:
    set file_id, "*.bam*" into indexed_bam_file

  script:
"""
sambamba index -t ${task.cpus} ${bam}
"""
}

bam_files_indexed = bam_files_bigwig.join(indexed_bam_file, by: 0)

process bam_to_bigwig {
  tag "$file_id"
  cpus 4
  publishDir "results/mapping/bigwig/", mode: 'copy'

  input:
    set file_id, file(bam), file(idx) from bam_files_indexed

  output:
    set file_id, "*.bw" into sorted_bam_files

  script:
"""
bamCoverage -p ${task.cpus} --ignoreDuplicates -b ${bam} -o ${file_id}.bw
"""
}

