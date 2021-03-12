version = "0.2.8"
container_url = "lbmc/g2gtools:${version}"

process vci_build {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.chain"), emit: vci
  script:
"""
g2gtools vcf2vci \
  -p ${task.cpus} \
  -f ${fasta} \
  -i ${vcf} \
  -s ${file_id} \
  -o ${file_id}.vci
"""
}