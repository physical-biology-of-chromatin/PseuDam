version = "2.2.1"
container_url = "lbmc/histat2:${version}"

params.index_fasta = ""
params.index_fasta_out = ""
process index_fasta {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.index_fasta_out != "") {
    publishDir "results/${params.index_fasta_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta)

  output:
    tuple val(file_id), path("*.ht2*"), emit: index
    tuple val(file_id), path("*_report.txt"), emit: report

  script:
"""
hisat2-build -p ${task.cpus} \
  ${fasta} \
  ${fasta.simpleName} &> \
  ${fasta.simpleName}_hisat2_index_report.txt

if grep -q "Error" ${fasta.simpleName}_bowtie2_index_report.txt; then
  exit 1
fi
"""
}

params.mapping_fastq = "--very-sensitive"
params.mapping_fastq_out = ""
process mapping_fastq {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.mapping_fastq_out != "") {
    publishDir "results/${params.mapping_fastq_out}", mode: 'copy'
  }

  input:
  tuple val(index_id), path(index)
  tuple val(file_id), path(reads)

  output:
  tuple val(file_id), path("*.bam"), emit: bam
  path "*_report.txt", emit: report

  script:
  index_id = index[0]
  for (index_file in index) {
    if (index_file =~ /.*\.1\.ht2.*/) {
        index_id = ( index_file =~ /(.*)\.1\.ht2.*/)[0][1]
    }
  }
  switch(file_id) {
    case {it instanceof List}:
      file_prefix = file_id[0]
    break
    case {it instanceof Map}:
      file_prefix = file_id.values()[0]
    break
    default:
      file_prefix = file_id
    break
  }

  if (reads.size() == 2)
  """
  histat2 ${params.mapping_fastq} \
    -p ${task.cpus} \
    -x ${index_id} \
    -1 ${reads[0]} \
    -2 ${reads[1]} 2> \
    ${file_prefix}_ht2_mapping_report_tmp.txt \
    | samtools view -@ ${task.cpus} -bS - \
    | samtools sort -@ ${task.cpus} -o ${file_prefix}.bam -

  if grep -q "Error" ${file_prefix}_ht2_mapping_report_tmp.txt; then
    exit 1
  fi
  """
  else
  """
  histat2 ${params.mapping_fastq} \
    -p ${task.cpus} \
    -x ${index_id} \
    -U ${reads} 2> \
    ${file_prefix}_ht2_mapping_report_tmp.txt \
    | samtools view -@ ${task.cpus} -bS - \
    | samtools sort -@ ${task.cpus} -o ${file_prefix}.bam -

  if grep -q "Error" ${file_prefix}_ht2_mapping_report_tmp.txt; then
    exit 1
  fi
  """
}
