params.fastq = "$baseDir/data/fastq/*.fastq"

log.info "fastq files : ${params.fastq}"
log.info "index files : ${params.index}"

Channel
  .fromPath( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { fastq_files }
Channel
  .fromPath( params.index )
  .ifEmpty { error "Cannot find any index files matching: ${params.index}" }
  .set { index_files }

process mapping_fastq {
  tag "$file_id"
  publishDir "results/mapping/bams/", mode: 'copy'

  input:
  set file_id, file(reads) from fastq_files
  file index from index_files.collect()

  output:
  set file_id, "*.bam" into bam_files
  file "*.out" into mapping_report

  script:
"""
mkdir -p index
mv ${index} index/
STAR --runThreadN ${task.cpus} \
--genomeDir index/ \
--readFilesIn ${reads} \
--outFileNamePrefix ${file_id} \
--outSAMmapqUnique 0 \
--outSAMtype BAM SortedByCoordinate
"""
}
