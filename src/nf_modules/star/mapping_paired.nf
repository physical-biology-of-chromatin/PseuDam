params.fastq = "$baseDir/data/fastq/*_{1,2}.fastq"
params.index = "$baseDir/data/index/*.index.*"

log.info "fastq files : ${params.fastq}"
log.info "index files : ${params.index}"

Channel
  .fromFilePairs( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }
Channel
  .fromPath( params.index )
  .ifEmpty { error "Cannot find any index files matching: ${params.index}" }
  .set { index_files }

process mapping_fastq {
  tag "$pair_id"
  publishDir "results/mapping/bams/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files
  file index from index_files.collect()

  output:
  set pair_id, "*.bam" into bam_files
  file "*.out" into mapping_report

  script:
"""
mkdir -p index
mv ${index} index/
STAR --runThreadN ${task.cpus} \
--genomeDir index/ \
--readFilesIn ${reads[0]} ${reads[1]} \
--outFileNamePrefix ${pair_id} \
--outSAMmapqUnique 0 \
--outSAMtype BAM SortedByCoordinate
"""
}

