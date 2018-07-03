params.fastq = "$baseDir/data/fastq/*_{1,2}.fastq"
params.fasta = "$baseDir/data/fasta/*.fasta"
params.bed = "$baseDir/data/annot/*.bed"

log.info "fasta file : ${params.fasta}"
log.info "bed file : ${params.bed}"
log.info "fastq files : ${params.fastq}"

Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
  .set { fasta_files }

Channel
  .fromPath( params.bed )
  .ifEmpty { error "Cannot find any bed files matching: ${params.bed}" }
  .set { bed_files }

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

process trimming {
  tag "${reads}"
  cpus 4
  publishDir "results/fastq/trimming/", mode: 'copy'

  input:
  file reads from fastq_files_cut

  output:
  file "*_trim_R{1,2}.fastq.gz" into fastq_files_trim

  script:
"""
UrQt --t 20 --m ${task.cpus} --gz \
--in ${reads[0]} --inpair ${reads[1]} \
--out ${reads[0].baseName}_trim_R1.fastq.gz --outpair ${reads[1].baseName}_trim_R2.fastq.gz \
> ${reads[0].baseName}_trimming_report.txt
"""
}

process fasta_from_bed {
  tag "${bed.baseName}"
  cpus 4
  publishDir "results/fasta/", mode: 'copy'

  input:
  file fasta from fasta_files
  file bed from bed_files

  output:
  file "*_extracted.fasta" into fasta_files_extracted

  script:
"""
bedtools getfasta -name \
-fi ${fasta} -bed ${bed} -fo ${bed.baseName}_extracted.fasta
"""
}

process index_fasta {
  tag "$fasta.baseName"
  publishDir "results/mapping/index/", mode: 'copy'

  input:
    file fasta from fasta_files_extracted

  output:
    file "*.index*" into index_files

  script:
"""
kallisto index -k 31 --make-unique -i ${fasta.baseName}.index ${fasta} \
> ${fasta.baseName}_kallisto_report.txt
"""
}


process mapping_fastq {
  tag "$reads"
  cpus 4
  publishDir "results/mapping/quantification/", mode: 'copy'

  input:
  file reads from fastq_files_trim
  file index from index_files

  output:
  file "*" into counts_files

  script:
"""
mkdir ${reads[0].baseName}
kallisto quant -i ${index} -t ${task.cpus} \
--bias --bootstrap-samples 100 -o ${reads[0].baseName} \
${reads[0]} ${reads[1]} &> ${reads[0].baseName}_kallisto_report.txt
"""
}




