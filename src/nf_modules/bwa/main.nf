version = "0.7.17"
container_url = "lbmc/bwa:${version}"

process index_fasta {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(fasta)

  output:
    tuple val(file_id), path("${fasta.simpleName}.*"), emit: index
    tuple val(file_id), path("*_bwa_report.txt"), emit: report

  script:
"""
bwa index -p ${fasta.simpleName} ${fasta} \
&> ${fasta.simpleName}_bwa_report.txt
"""
}


process mapping_fastq {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
  tuple val(file_id), path(reads)
  tuple val(index_id), path(index)

  output:
  tuple val(file_id), path("*.bam"), emit: bam
  tuple val(file_id), path("${file_id}_bwa_report.txt"), emit: report

  script:
if (file_id.containsKey('library')) {
  library = file_id.library
} else {
  library = file_id
}
bwa_mem_R = "@RG\\tID:${library}\\tSM:${library}\\tLB:lib_${library}\\tPL:illumina"
"""
bwa mem -t ${task.cpus} \
-R '${params.bwa_mem_R}' \
${index_id} ${reads[0]} ${reads[1]} 2> \
  ${file_id}_bwa_report.txt | \
  samtools view -@ ${task.cpus} -Sb - > ${file_id}.bam
"""
}

