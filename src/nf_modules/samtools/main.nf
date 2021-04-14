version = "1.11"
container_url = "lbmc/samtools:${version}"

process index_fasta {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(fasta)
  output:
    tuple val(file_id), path("*.fai"), emit: index

  script:
"""
samtools faidx ${fasta}
"""
}

filter_bam_quality_threshold = 30

process filter_bam_quality {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*_filtered.bam"), emit: bam
  script:
"""
samtools view -@ ${task.cpus} -hb ${bam} -q ${filter_bam_quality_threshold} > \
  ${bam.simpleName}_filtered.bam
"""
}


process filter_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam)
    tuple val(bed_id), path(bed)

  output:
    tuple val(file_id), path("*_filtered.bam"), emit: bam
  script:
"""
samtools view -@ ${task.cpus} -hb ${bam} -L ${bed} > \
  ${bam.simpleName}_filtered.bam
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
samtools view -@ ${task.cpus} -F 4 -hb ${bam} > \
  ${bam.simpleName}_mapped.bam
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
samtools view -@ ${task.cpus} -f 4 -hb ${bam} > ${bam.simpleName}_unmapped.bam
"""
}


process index_bam {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path(bam), emit: bam
    tuple val(file_id), path("*.bam.bai"), emit: bam_idx

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
samtools sort -@ ${task.cpus} -O BAM -o ${bam.simpleName}_sorted.bam ${bam}
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
samtools view --@ ${Math.round(task.cpus/2)} \
  -hb -F 0x10 ${bam} > ${bam.simpleName}_forward.bam &
samtools view --@ ${Math.round(task.cpus/2)} \
  -hb -f 0x10 ${bam} > ${bam.simpleName}_reverse.bam
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
samtools merge ${first_bam} ${second_bam} \
  ${first_bam.simpleName}_${second_file.simpleName}.bam
"""
}

process merge_multi_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  cpus = 2

  input:
    tuple val(file_id), path(bams)

  output:
    tuple val(file_id), path("*_merged.bam*"), emit: bam
  script:
"""
samtools merge -@ ${task.cpus} \
  ${bams[0].simpleName}_merged.bam \
  ${bams}
"""
}

process stats_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  cpus = 2

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*.tsv"), emit: tsv
  script:
"""
samtools flagstat -@ ${task.cpus} -O tsv ${bam} > ${bam.simpleName}_stats.tsv
"""
}

process flagstat_2_multiqc {
  tag "$file_id"

  input:
    tuple val(file_id), path(tsv)

  output:
    tuple val(file_id), path("*.txt"), emit: report
"""
mv ${tsv} ${tsv.simpleName}.flagstat.txt
"""
}

process idxstat_2_multiqc {
  tag "$file_id"

  input:
    tuple val(file_id), path(tsv)

  output:
    tuple val(file_id), path("*.txt"), emit: report
"""
mv ${tsv} ${tsv.simpleName}.idxstats.txt
"""
}