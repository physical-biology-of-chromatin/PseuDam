version = "0.6.7"
container_url = "lbmc/sambamba:${version}"

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

process sort_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*.bam*"), emit: bam

  script:
"""
sambamba sort -t ${task.cpus} -o ${file_id}_sorted.bam ${bam}
"""
}


process split_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*_forward.bam*"), emit: bam_forward
    tuple val(file_id), path("*_reverse.bam*"), emit: bam_reverse
  script:
"""
sambamba view -t ${task.cpus} -h -F "strand == '+'" ${bam} > ${file_id}_forward.bam
sambamba view -t ${task.cpus} -h -F "strand == '-'" ${bam} > ${file_id}_reverse.bam
"""
}
