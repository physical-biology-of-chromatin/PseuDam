params.fastq = "$baseDir/data/fastq/*.fastq"

log.info "fastq files : ${params.fastq}"
log.info "index files : ${params.index}"

Channel
  .fromPath( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }
Channel
  .fromPath( params.index )
  .ifEmpty { error "Cannot find any index files matching: ${params.index}" }
  .set { index_files }

process mapping_fastq {
  tag "$reads.baseName"
  cpus 4
  publishDir "results/mapping/bams/", mode: 'copy'

  input:
  file reads from fastq_files
  file index from index_files.toList()

  output:
  file "*.bam" into bam_files
  file "*_report.txt" into mapping_report

  script:
index_id = index[0]
for (index_file in index) {
  if (index_file =~ /.*\.1\.ebwt/ && !(index_file =~ /.*\.rev\.1\.ebwt/)) {
      index_id = ( index_file =~ /(.*)\.1\.ebwt/)[0][1]
  }
}
"""
bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_id} \
-q ${reads} 2> \
${reads.baseName}_bowtie_report.txt | \
samtools view -Sb - > ${reads.baseName}.bam

if grep -q "Error" ${reads.baseName}_bowtie_report.txt; then
  exit 1
fi
"""
}
