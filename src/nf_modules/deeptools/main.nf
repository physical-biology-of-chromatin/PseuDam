version = "3.0.2"
container_url = "lbmc/deeptools:${version}"

process index_bam {
  container = "${container_url}"
  label "big_mem__cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*.bam*"), emit: bam

  script:
"""
sambamba index -t ${task.cpus} ${bam}
"""
}

process bam_to_bigwig {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  publishDir "results/mapping/bigwig/", mode: 'copy'

  input:
    tuple val(file_id), path(bam), path(idx)

  output:
    tuple val(file_id), path("*.bw"), emit: bw

  script:
"""
bamCoverage -p ${task.cpus} --ignoreDuplicates -b ${bam} -o ${file_id}.bw
"""
}

process compute_matrix {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "${bed_file_id}"
  publishDir "results/metagene/region_matrix/", mode: 'copy'

  input:
    tuple val(file_id), path(bw)
    tuple val(bed_file_id), path(bed)

  output:
    tuple val(bed_file_id), path("*.mat.gz"), emit: matrix

  script:
"""
computeMatrix scale-regions -S ${bw} \
  -p ${task.cpus} \
  -R ${bed} \
  --beforeRegionStartLength 100 \
  --afterRegionStartLength 100 \
  -o ${bed_file_id}.mat.gz
"""
}

process plot_profile {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  publishDir "results/metagene/region_matrix/", mode: 'copy'

  input:
    tuple val(file_id), path(matrix)

  output:
    tuple val(file_id), path("*.pdf"), emit: pdf

  script:
/*
see more option at
https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html
*/
"""
plotProfile -m ${matrix} \
  --plotFileFormat=pdf \
  -out ${file_id}.pdf \
  --plotType=fill \
  --perGroup \
  --plotTitle "${params.title}"
"""
}
