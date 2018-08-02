params.fastq = "$baseDir/data/fastq/*.fastq"
params.index = "$baseDir/data/index/*.index*"
params.mean = 125
params.sd = 100

log.info "fastq files : ${params.fastq}"
log.info "index files : ${params.index}"
log.info "mean read size: ${params.mean}"
log.info "sd read size: ${params.sd}"

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
  publishDir "results/mapping/quantification/", mode: 'copy'

  input:
  file reads from fastq_files
  file index from index_files.collect()

  output:
  file "*" into count_files

  script:
index_name = (index[0].baseName =~ /(.*)\.\d/)[0][1]
"""
rsem-calculate-expression --bowtie2 \
--bowtie2-path \$(which bowtie2 | sed 's/bowtie2\$//g') \
--bowtie2-sensitivity-level "very_sensitive" \
--fragment-length-mean ${params.mean} --fragment-length-sd ${params.sd} \
--output-genome-bam -p ${task.cpus} \
${reads} ${index_name} ${reads.baseName} \
> ${reads.baseName}_rsem_bowtie2_report.txt
"""
}

