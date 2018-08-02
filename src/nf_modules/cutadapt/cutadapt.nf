/*
* cutadapt :
* Imputs : fastq files
* Output : fastq files
*/

/*                      Illumina adaptor removal                             */

/*
* for paired-end data
*/

params.fastq = "$baseDir/data/fastq/*_{1,2}.fastq"

log.info "fastq files : ${params.fastq}"

Channel
  .fromFilePairs( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

process adaptor_removal {
  tag "$pair_id"
  publishDir "results/fastq/adaptor_removal/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files

  output:
  set pair_id, "*_cut_R{1,2}.fastq.gz" into fastq_files_cut

  script:
  """
  cutadapt -a AGATCGGAAGAG -g CTCTTCCGATCT -A AGATCGGAAGAG -G CTCTTCCGATCT \
  -o ${pair_id}_cut_R1.fastq.gz -p ${pair_id}_cut_R2.fastq.gz \
  ${reads[0]} ${reads[1]} > ${pair_id}_report.txt
  """
}

/*
* for single-end data
*/

log.info "fastq files : ${params.fastq}"

Channel
  .fromPath( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { fastq_files }

process adaptor_removal {
  tag "$file_id"

  input:
  set file_id, file(reads) from fastq_files

  output:
  set file_id, "*_cut.fastq.gz" into fastq_files_cut

  script:
  """
  cutadapt -a AGATCGGAAGAG -g CTCTTCCGATCT\
  -o ${file_id}_cut.fastq.gz \
  ${reads} > ${file_id}_report.txt
  """
}


/*                      quality trimming                                     */

/*
* for paired-end data
*/

params.fastq = "$baseDir/data/fastq/*_{1,2}.fastq"

log.info "fastq files : ${params.fastq}"

Channel
  .fromFilePairs( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

process trimming {
  tag "$pair_id"
  publishDir "results/fastq/trimming/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files

  output:
  set pair_id, "*_trim_R{1,2}.fastq.gz" into fastq_files_trim

  script:
  """
  cutadapt -q 20,20 \
  -o ${pair_id}_trim_R1.fastq.gz -p ${pair_id}_trim_R2.fastq.gz \
  ${reads[0]} ${reads[1]} > ${pair_id}_report.txt
  """
}

/*
* for single-end data
*/

log.info "fastq files : ${params.fastq}"

Channel
  .fromPath( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { fastq_files }

process trimming {
  tag "$file_id"

  input:
  set file_id, file(reads) from fastq_files

  output:
  set file_id, "*_trim.fastq.gz" into fastq_files_cut

  script:
  """
  cutadapt -q 20,20 \
  -o ${file_id}_trim.fastq.gz \
  ${reads} > ${file_id}_report.txt
  """
}

