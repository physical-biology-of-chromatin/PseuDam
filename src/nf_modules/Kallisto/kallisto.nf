/*
* Kallisto :
* Imputs : fastq files
* Imputs : fasta files
* Output : bam files
*/

/*                      fasta indexing                                     */
params.fasta = "$baseDir/data/bam/*.fasta"

log.info "fasta files : ${params.fasta}"

Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
  .set { fasta_file }

process index_fasta {
  tag "$fasta.baseName"
  publishDir "results/mapping/index/", mode: 'copy'

  input:
    file fasta from fasta_file

  output:
    file "*.index*" into index_files

  script:
"""
kallisto index -k 31 --make-unique -i ${fasta.baseName}.index ${fasta} \
> ${fasta.baseName}_kallisto_report.txt
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
  tag "$reads"
  cpus 4
  publishDir "results/mapping/quantification/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files
  file index from index_files.collect()

  output:
  file "*" into counts_files

  script:
"""
mkdir ${reads[0].baseName}
kallisto quant -i ${index} -t ${task.cpus} \
--bias --bootstrap-samples 100 -o ${pair_id} \
${reads[0]} ${reads[1]} &> ${pair_id}_kallisto_report.txt
"""
}


/*
* for single-end data
*/

params.fastq = "$baseDir/data/fastq/*.fastq"
params.index = "$baseDir/data/index/*.index*"
params.mean = 200
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
"""
mkdir ${file_id}
kallisto quant -i ${index} -t ${task.cpus} --single \
--bias --bootstrap-samples 100 -o ${file_id} \
-l ${params.mean} -s ${params.sd} \
${reads} > ${file_id}_kallisto_report.txt
"""
}


