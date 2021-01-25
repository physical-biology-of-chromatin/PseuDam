version = "1.7"
container_url = "lbmc/samtools:${version}"

process filter_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam)
    path bed

  output:
    tuple val(file_id), path("*_filtered.bam"), emit: bam
  script:
"""
samtools view -@ ${task.cpus} -hb ${bam} -L ${bed} > ${file_id}_filtered.bam
"""
}

process filter_bam_mapped {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*_mapped.bam"), emit: bam
  script:
"""
samtools view -@ ${task.cpus} -F 4 -hb ${bam} -L ${bed} > ${file_id}_mapped.bam
"""
}

process filter_bam_unmapped {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*_unmapped.bam"), emit: bam
  script:
"""
samtools view -@ ${task.cpus} -f 4 -hb ${bam} > ${file_id}_unmapped.bam
"""
}


process index_bam {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*.bam*"), emit: bam

  script:
"""
samtools index ${bam}
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
samtools sort -@ ${task.cpus} -O BAM -o ${file_id}_sorted.bam ${bam}
"""
}


process split_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  cpus = 2

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*_forward.bam*"), emit: bam_forward
    tuple val(file_id), path("*_reverse.bam*"), emit: bam_reverse
  script:
"""
samtools view -hb -F 0x10 ${bam} > ${file_id}_forward.bam &
samtools view -hb -f 0x10 ${bam} > ${file_id}_reverse.bam
"""
}


process merge_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  cpus = 2

  input:
    tuple val(first_file_id), path(first_bam)
    tuple val(second_file_id), path(second_bam)

  output:
    tuple val(file_id), path("*.bam*"), emit: bam
  script:
"""
samtools merge ${first_bam} ${second_bam} ${first_bam_id}_${second_file_id}.bam
"""
}
