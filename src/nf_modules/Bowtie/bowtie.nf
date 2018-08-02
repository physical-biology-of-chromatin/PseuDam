/*
* Bowtie :
* Imputs : fastq files
* Imputs : fasta files
* Output : bam files
*/

/*                      fasta indexing                                     */
params.fasta = "$baseDir/data/bam/*.fasta"

log.info "fasta files : ${params.fasta}"

Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any bam files matching: ${params.fasta}" }
  .set { fasta_file }

process index_fasta {
  tag "$fasta.baseName"
  cpus 4
  publishDir "results/mapping/index/", mode: 'copy'

  input:
    file fasta from fasta_file

  output:
    file "*.index*" into index_files
    file "*_report.txt" into indexing_report

  script:
"""
bowtie-build --threads ${task.cpus} -f ${fasta} ${fasta.baseName}.index &> ${fasta.baseName}_bowtie_report.txt

if grep -q "Error" ${fasta.baseName}_bowtie_report.txt; then
  exit 1
fi
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
  publishDir "results/mapping/bams/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files
  file index from index_files.collect()

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
# -v specify the max number of missmatch, -k the number of match reported per
# reads
bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_id} \
-1 ${reads[0]} -2 ${reads[1]} 2> \
${pair_id}_bowtie_report.txt | \
samtools view -Sb - > ${pair_id}.bam

if grep -q "Error" ${pair_id}_bowtie_report.txt; then
  exit 1
fi
"""
}


/*
* for single-end data
*/
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

