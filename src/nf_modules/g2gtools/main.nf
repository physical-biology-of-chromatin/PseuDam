version = "0.2.8"
container_url = "lbmc/g2gtools:${version}"

process vci_build {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta)
  output:
    tuple val(file_id), path("*.vci.gz"), emit: vci
  script:
"""
g2gtools vcf2vci \
  -p ${task.cpus} \
  -f ${fasta} \
  -i ${vcf} \
  -s ${file_id} \
  -o ${file_id}.vci.gz
"""
}

process incorporate_snp {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(vci)
    tuple val(ref_id), path(fasta)
  output:
    tuple val(file_id), path("${file_id}_snp.fasta"), path("${vci}"), emit: fasta
  script:
"""
g2gtools patch \
  -p ${task.cpus} \
  -i ${fasta} \
  -c ${vci} \
  -o ${file_id}_snp.fasta
"""
}

process incorporate_indel {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(fasta), path(vci)
  output:
    tuple val(file_id), path("${file_id}_snp_indel.fasta"), path("${vci}"), emit: fasta
  script:
"""
g2gtools transform \
  -p ${task.cpus} \
  -i ${fasta} \
  -c ${vci} \
  -o ${file_id}_snp_indel.fasta
"""
}