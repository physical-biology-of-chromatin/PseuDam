params.fastq = "$baseDir/data/*.fastq"
params.fasta = "$baseDir/data/*.fasta"
log.info "fastq files : ${params.fastq}"
log.info "fasta files : ${params.fasta}"

Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any bam files matching: ${params.fasta}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { fasta_file }
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

process trimming {
  tag "${reads}"
  cpus 4
  publishDir "results/fastq/trimming/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files_cut

  output:
  set pair_id, "*_trim_R{1,2}.fastq.gz" into fastq_files_trim

  script:
"""
UrQt --t 20 --m ${task.cpus} --gz \
--in ${reads[0]} --inpair ${reads[1]} \
--out ${pair_id}_trim_R1.fastq.gz --outpair ${pair_id}_trim_R2.fastq.gz \
> ${pair_id}_trimming_report.txt
"""
}

process index_fasta {
  tag "$fasta_id"
  cpus 4
  publishDir "results/mapping/index/", mode: 'copy'

  input:
    set fasta_id, file(fasta) from fasta_file

  output:
    set fasta_id, "${fasta.baseName}.*" into index_files
    file "*_bwa_report.txt" into index_files_report

  script:
"""
bwa index -p ${fasta.baseName} ${fasta} \
&> ${fasta.baseName}_bwa_report.txt
"""
}


process mapping_fastq {
  tag "$reads"
  cpus 4
  publishDir "results/mapping/sam/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files_trim
  set index_id, file(index) from index_files.collect()

  output:
  file "${pair_id}.sam" into sam_files
  file "${pair_id}_bwa_report.txt" into mapping_repport_files

  script:
"""
bwa mem -t ${task.cpus} \
${index_id} ${reads[0]} ${reads[1]} \
-o ${pair_id}.sam &> ${pair_id}_bwa_report.txt
"""
}

