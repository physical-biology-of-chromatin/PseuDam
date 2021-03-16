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
  -s ${file_id.library} \
  -o ${file_id.id}.vci.gz
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
    tuple val(file_id), path("${file_id.id}_snp.fasta"), path("${vci}"), emit: fasta
  script:
"""
g2gtools patch \
  -p ${task.cpus} \
  -i ${fasta} \
  -c ${vci} \
  -o ${file_id.id}_snp.fasta
"""
}

process incorporate_indel {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(fasta), path(vci)
  output:
    tuple val(file_id), path("${file_id.id}_snp_indel.fasta"), path("${vci}"), emit: fasta
  script:
"""
g2gtools transform \
  -p ${task.cpus} \
  -i ${fasta} \
  -c ${vci} \
  -o ${file_id.id}_snp_indel.fasta
"""
}

process convert_gtf {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(vci)
    tuple val(annot_id), path(gtf)
  output:
    tuple val(file_id), path("${file_id.id}.gtf"), emit: gtf
  script:
"""
g2gtools convert \
  -i ${gtf} \
  -c ${vci} \
  -o ${file_id.id}.gtf
"""
}

process convert_bed {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(vci)
    tuple val(annot_id), path(bed)
  output:
    tuple val(file_id), path("${file_id.id}.bed"), emit: bed
  script:
"""
g2gtools convert \
  -i ${bed} \
  -c ${vci} \
  -o ${file_id.id}.bed
"""
}

process convert_bam {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${bam_id} ${file_id}"

  input:
    tuple val(file_id), path(vci)
    tuple val(bam_id), path(bam)
  output:
    tuple val(file_id), path("${file_id.id}_${bam_id.baseName}.bam"), emit: bam
  script:
"""
g2gtools convert \
  -i ${bam} \
  -c ${vci} \
  -o ${file_id.id}_${bam.baseName}.bam
"""
}