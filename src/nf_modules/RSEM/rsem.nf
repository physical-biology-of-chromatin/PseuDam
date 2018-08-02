/*
* RSEM :
* Imputs : fastq files
* Imputs : fasta files
* Output : bam files
*/

/*                      fasta indexing                                     */
params.fasta = "$baseDir/data/bam/*.fasta"
params.annotation = "$baseDir/data/bam/*.gff3"

log.info "fasta files : ${params.fasta}"

Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
  .set { fasta_file }
Channel
  .fromPath( params.annotation )
  .ifEmpty { error "Cannot find any annotation files matching: ${params.annotation}" }
  .set { annotation_file }

process index_fasta {
  tag "$fasta.baseName"
  cpus 4
  publishDir "results/mapping/index/", mode: 'copy'

  input:
    file fasta from fasta_file
    file annotation from annotation_file

  output:
    file "*.index*" into index_files

  script:
  def cmd_annotation = "--gff3 ${annotation}"
  if(annotation ==~ /.*\.gtf$/){
    cmd_annotation = "--gtf ${annotation}"
  }
"""
rsem-prepare-reference -p ${task.cpus} --bowtie2 \
--bowtie2-path \$(which bowtie2 | sed 's/bowtie2\$//g') \
${cmd_annotation} ${fasta} ${fasta.baseName}.index > \
${fasta.baseName}_rsem_bowtie2_report.txt
"""
}


/*
* for paired-end data
*/
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
  cpus 4
  publishDir "results/mapping/quantification/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files
  file index from index_files.collect()

  output:
  file "*" into counts_files

  script:
index_name = (index[0].baseName =~ /(.*)\.\d/)[0][1]
"""
rsem-calculate-expression --bowtie2 \
--bowtie2-path \$(which bowtie2 | sed 's/bowtie2\$//g') \
--bowtie2-sensitivity-level "very_sensitive" \
-output-genome-bam -p ${task.cpus} \
--paired-end ${reads[0]} ${reads[1]} ${index_name} ${pair_id} \
> ${pair_id}_rsem_bowtie2_report.txt
"""
}



/*
* for single-end data
*/

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
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { fastq_files }
Channel
  .fromPath( params.index )
  .ifEmpty { error "Cannot find any index files matching: ${params.index}" }
  .set { index_files }

process mapping_fastq {
  tag "$file_id"
  cpus 4
  publishDir "results/mapping/quantification/", mode: 'copy'

  input:
  set file_id, file(reads) from fastq_files
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
${reads} ${index_name} ${file_id} \
> ${reads.baseName}_rsem_bowtie2_report.txt
"""
}
