version = "2.17"
container_url = "lbmc/minimap2:${version}"

params.index_fasta = ""
process index_fasta {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$fasta.baseName"

  input:
    path fasta

  output:
    tuple path("${fasta}"), path("*.mmi*"), emit: index
    path "*_report.txt", emit: report

  script:
  memory = "${task.memory}" - ~/\s*GB/
"""
minimap2 ${params.index_fasta} -t ${task.cpus} -I ${memory}G -d ${fasta.baseName}.mmi ${fasta}
"""
}

params.mapping_fastq = "-ax sr"
process mapping_fastq {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$pair_id"

  input:
  tuple path(fasta), path(index)
  tuple val(pair_id), path(reads)

  output:
  tuple val(pair_id), path("*.bam"), emit: bam
  path "*_report.txt", emit: report

  script:
  memory = "${task.memory}" - ~/\s*GB/
  memory = memory / (task.cpus + 1.0)
if (reads instanceof List)
"""
minimap2 ${params.mapping_fastq} -t ${task.cpus} -K ${memory} ${fasta} ${reads[0]} ${reads[1]} |
  samtools view -Sb - > ${pair_id}.bam
"""
else
"""
minimap2 ${params.mapping_fastq} -t ${task.cpus} -K ${memory} ${fasta} ${reads} |
  samtools view -Sb - > ${reads.baseName}.bam
"""
}