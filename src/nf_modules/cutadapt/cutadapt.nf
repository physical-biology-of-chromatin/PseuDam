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
  file "*_cut_R{1,2}.fastq.gz" into fastq_files_cut

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

params.fastq = "$baseDir/data/fastq/*.fastq"

log.info "fastq files : ${params.fastq}"

Channel
  .fromPath( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

process adaptor_removal {
  tag "$reads.baseName"
  publishDir "results/fastq/adaptor_removal/", mode: 'copy'

  input:
  file reads from fastq_files

  output:
  file "*_cut.fastq.gz" into fastq_files_cut

  script:
  """
  cutadapt -a AGATCGGAAGAG -g CTCTTCCGATCT\
  -o ${reads.baseName}_cut.fastq.gz \
  ${reads} > ${reads.baseName}_report.txt
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
  file "*_trim_R{1,2}.fastq.gz" into fastq_files_trim

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

params.fastq = "$baseDir/data/fastq/*.fastq"

log.info "fastq files : ${params.fastq}"

Channel
  .fromPath( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

process trimming {
  tag "$reads.baseName"
  publishDir "results/fastq/trimming/", mode: 'copy'

  input:
  file reads from fastq_files

  output:
  file "*_trim.fastq.gz" into fastq_files_trim

  script:
  """
  cutadapt -q 20,20 \
  -o ${reads.baseName}_trim.fastq.gz \
  ${reads} > ${reads.baseName}_report.txt
  """
}

