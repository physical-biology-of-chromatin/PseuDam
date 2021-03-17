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
    tuple val(file_id), path("*_report.txt"), emit: report
  script:
  input_vcf = ""
  for (vcf_file in vcf) {
    input_vcf += " -i ${vcf_file}"
  }
"""
g2gtools vcf2vci \
  -p ${task.cpus} \
  -f ${fasta} \
  ${input_vcf} \
  -s ${file_id} \
  -o ${file_id}.vci.gz 2> ${file_id}_g2gtools_vcf2vci_report.txt
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
    tuple val(file_id), path("*_report.txt"), emit: report
  script:
"""
g2gtools patch \
  -p ${task.cpus} \
  -i ${fasta} \
  -c ${vci} \
  -o ${file_id}_snp.fasta 2> ${file_id}_g2gtools_path_report.txt
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
    tuple val(file_id), path("*_report.txt"), emit: report
  script:
"""
g2gtools transform \
  -p ${task.cpus} \
  -i ${fasta} \
  -c ${vci} \
  -o ${file_id}_snp_indel.fasta 2> ${file_id}_g2gtools_transform_report.txt
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
    tuple val(file_id), path("${file_id}.gtf"), emit: gtf
    tuple val(file_id), path("*_report.txt"), emit: report
  script:
"""
g2gtools convert \
  -i ${gtf} \
  -c ${vci} \
  -o ${file_id}.gtf 2> ${file_id}_g2gtools_convert_report.txt
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
    tuple val(file_id), path("${file_id}.bed"), emit: bed
    tuple val(file_id), path("*_report.txt"), emit: report
  script:
"""
g2gtools convert \
  -i ${bed} \
  -c ${vci} \
  -o ${file_id}.bed 2> ${file_id}_g2gtools_convert_report.txt
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
    tuple val(file_id), path("${file_id}_${bam_id.baseName}.bam"), emit: bam
    tuple val(file_id), path("*_report.txt"), emit: report
  script:
"""
g2gtools convert \
  -i ${bam} \
  -c ${vci} \
  -o ${file_id}_${bam.baseName}.bam 2> ${file_id}_g2gtools_convert_report.txt
"""
}