version = "0.44.0"
container_url = "lbmc/kallisto:${version}"

process index_fasta {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$fasta.baseName"
  publishDir "results/mapping/index/", mode: 'copy'

  input:
    path fasta

  output:
    path "*.index*", emit: index
    path "*_report.txt", emit: report

  script:
"""
kallisto index -k 31 --make-unique -i ${fasta.baseName}.index ${fasta} \
2> ${fasta.baseName}_kallisto_index_report.txt
"""
}


process mapping_fastq_pairedend {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$pair_id"
  publishDir "results/mapping/counts/", mode: 'copy'

  input:
  path index
  tuple val(pair_id), path(reads)

  output:
  path "${pair_id}", emit: counts
  path "*_report.txt", emit: report

  script:
"""
mkdir ${pair_id}
kallisto quant -i ${index} -t ${task.cpus} \
--bias --bootstrap-samples 100 -o ${pair_id} \
${reads[0]} ${reads[1]} &> ${pair_id}_kallisto_mapping_report.txt
"""
}


process mapping_fastq_singleend {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  publishDir "results/mapping/counts/", mode: 'copy'

  input:
  path index
  tuple val(file_id), path(reads)

  output:
  tuple val(file_id), path("${pair_id}"), emit: counts
  path "*_report.txt", emit: report

  script:
"""
mkdir ${file_id}
kallisto quant -i ${index} -t ${task.cpus} --single \
--bias --bootstrap-samples 100 -o ${file_id} \
-l ${params.mean} -s ${params.sd} \
${reads} &> ${reads.simpleName}_kallisto_mapping_report.txt
"""
}
