version = "v2.2.2_cv3"
container_url = "biocontainers/danpos:${version}"

include {
  bigwig2_to_wig2
} from "./../ucsc/main.nf"

params.dpos = "--smooth_width 0 -n N "
params.dpos_out = ""

process dpos_bam {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.dpos_out != "") {
    publishDir "results/${params.dpos_out}", mode: 'copy'
  }

  input:
    tuple val(fastq_id), path(fastq)
    tuple val(file_id), path(bam_ip), path(bam_wce)

  output:
    tuple val(file_id), path("${file_prefix}/pooled/*.wig"), emit: wig
    tuple val(file_id), path("${file_prefix}"), emit: folder

  script:

  switch(file_id) {
    case {it instanceof List}:
      file_prefix = file_id[0]
    break
    case {it instanceof Map}:
      file_prefix = file_id.values()[0]
    break
    default:
      file_prefix = file_id
    break
  }

  m = 0
  if (fastq.size() == 2){
    m = 1
  }
"""
danpos.py dpos -m ${m} \
  ${params.dpos} \
  -b ${bam_wce} \
  -o ${file_prefix} \
  ${bam_ip}
"""
}

workflow dpos_bw {
  take:
    fastq
    bw
  main:
    dpos_wig(fastq, bigwig2_to_wig2(bw))
  emit:
  wig = dpos_wig.out.wig
  folder = dpos_wig.out.folder
}

process dpos_wig {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.dpos_out != "") {
    publishDir "results/${params.dpos_out}", mode: 'copy'
  }

  input:
    tuple val(fastq_id), path(fastq)
    tuple val(file_id), path(wig_ip), path(wig_wce)

  output:
    tuple val(file_id), path("${file_prefix}/pooled/*.wig"), emit: wig
    tuple val(file_id), path("${file_prefix}"), emit: folder

  script:

  switch(file_id) {
    case {it instanceof List}:
      file_prefix = file_id[0]
    break
    case {it instanceof Map}:
      file_prefix = file_id.values()[0]
    break
    default:
      file_prefix = file_id
    break
  }

  m = 0
  if (fastq.size() == 2){
    m = 1
  }
"""
danpos.py dpos -m ${m} \
  ${params.dpos} \
  -b ${wig_wce} \
  -o ${file_prefix} \
  ${wig_ip}
"""
}

workflow dwig_bwvsbw {
  take:
    fastq
    bw_a
    bw_b
  main:
    dpos_wigvswig(
      fastq,
      bigwig2_to_wig2(bw_a),
      bigwig2_to_wig2(bw_b),
    )
  emit:
  wig = dpeak_wig.out.wig
  folder = dpeak_wig.out.folder
}

process dpos_wigvswig {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.dpos_out != "") {
    publishDir "results/${params.dpos_out}", mode: 'copy'
  }

  input:
    tuple val(fastq_id), path(fastq)
    tuple val(file_id_a), path(wig_ip_a), path(wig_wce_a)
    tuple val(file_id_b), path(wig_ip_b), path(wig_wce_b)

  output:
    tuple val(file_id), path("${file_prefix}/pooled/*.wig"), emit: wig
    tuple val(file_id_a), path("${file_prefix}"), emit: folder

  script:

  switch(file_id_a) {
    case {it instanceof List}:
      file_prefix = file_id_a[0]
    break
    case {it instanceof Map}:
      file_prefix = file_id_a.values()[0]
    break
    default:
      file_prefix = file_id_a
    break
  }

  m = 0
  if (fastq.size() == 2){
    m = 1
  }
"""
danpos.py dpos -m ${m} \
  ${params.dpos} \
  -b ${wig_ip_a}:${wig_wce_a},${wig_ip_b}:${wig_wce_b} \
  -o ${file_prefix} \
  ${wig_ip_a}:${wig_ip_b}
"""
}

params.dpeak = "--smooth_width 0 -n N "
params.dpeak_out = ""

process dpeak_bam {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.dpeak_out != "") {
    publishDir "results/${params.dpeak_out}", mode: 'copy'
  }

  input:
    tuple val(fastq_id), path(fastq)
    tuple val(file_id), path(bam_ip), path(bam_wce)

  output:
    tuple val(file_id), path("${file_prefix}/*.wig"), emit: wig
    tuple val(file_id), path("${file_prefix}/*.bed"), emit: bed
    tuple val(file_id), path("${file_prefix}"), emit: folder

  script:

  switch(file_id) {
    case {it instanceof List}:
      file_prefix = file_id[0]
    break
    case {it instanceof Map}:
      file_prefix = file_id.values()[0]
    break
    default:
      file_prefix = file_id
    break
  }

  m = 0
  if (fastq.size() == 2){
    m = 1
  }
"""
danpos.py dpeak -m ${m} \
  ${params.dpeak} \
  -b ${bam_wce} \
  -o ${file_prefix} \
  ${bam_ip}
mv ${file_prefix}/pooled/* ${file_prefix}/
rm -R ${file_prefix}/pooled
awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2-1, \$3, "Interval_"NR-1, \$6, "+" }' ${file_prefix}/${bam_ip.simpleName}.bgsub.peaks.xls > ${file_prefix}/${bam_ip.simpleName}.bgsub.positions.bed
awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$4-1, \$4, "Interval_"NR-1, \$6, "+" }' ${file_prefix}/${bam_ip.simpleName}.bgsub.peaks.xls > ${file_prefix}/${bam_ip.simpleName}.bgsub.positions.summit.bed
"""
}

workflow dpeak_bw {
  take:
    fastq
    bw
  main:
    dpeak_wig(fastq, bigwig2_to_wig2(bw))
  emit:
  wig = dpeak_wig.out.wig
  folder = dpeak_wig.out.folder
}


process dpeak_wig {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.dpeak_out != "") {
    publishDir "results/${params.dpeak_out}", mode: 'copy'
  }

  input:
    tuple val(fastq_id), path(fastq)
    tuple val(file_id), path(wig_ip), path(wig_wce)

  output:
  tuple val(file_id), path("${file_prefix}/*.wig"), emit: wig
  tuple val(file_id), path("${file_prefix}/*.bed"), emit: bed
  tuple val(file_id), path("${file_prefix}"), emit: folder

  script:

  switch(file_id) {
    case {it instanceof List}:
      file_prefix = file_id[0]
    break
    case {it instanceof Map}:
      file_prefix = file_id.values()[0]
    break
    default:
      file_prefix = file_id
    break
  }

  m = 0
  if (fastq.size() == 2){
    m = 1
  }
"""
danpos.py dpeak -m ${m} \
  ${params.dpeak} \
  -b ${wig_wce} \
  -o ${file_prefix} \
  ${wig_ip}
mv ${file_prefix}/pooled/* ${file_prefix}/
rm -R ${file_prefix}/pooled
awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2-1, \$3, "Interval_"NR-1, \$6, "+" }' ${file_prefix}/${wig_ip.simpleName}.bgsub.peaks.xls > ${file_prefix}/${wig_ip.simpleName}.bgsub.positions.bed
awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$4-1, \$4, "Interval_"NR-1, \$6, "+" }' ${file_prefix}/${wig_ip.simpleName}.bgsub.peaks.xls > ${file_prefix}/${wig_ip.simpleName}.bgsub.positions.summit.bed
"""
}

workflow dpeak_bwvsbw {
  take:
    fastq
    bw_a
    bw_b
  main:
    dpeak_wigvswig(
      fastq,
      bigwig2_to_wig2(bw_a),
      bigwig2_to_wig2(bw_b),
    )
  emit:
  wig = dpeak_wig.out.wig
  folder = dpeak_wig.out.folder
}


process dpeak_wigvswig {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.dpeak_out != "") {
    publishDir "results/${params.dpeak_out}", mode: 'copy'
  }

  input:
    tuple val(fastq_id), path(fastq)
    tuple val(file_id_a), path(wig_ip_a), path(wig_wce_a)
    tuple val(file_id_b), path(wig_ip_b), path(wig_wce_b)

  output:
  tuple val(file_id_a), path("${file_prefix}/pooled/*.wig"), emit: wig
  tuple val(file_id_a), path("${file_prefix}"), emit: folder

  script:

  switch(file_id_a) {
    case {it instanceof List}:
      file_prefix = file_id_a[0]
    break
    case {it instanceof Map}:
      file_prefix = file_id_a.values()[0]
    break
    default:
      file_prefix = file_id_a
    break
  }

  m = 0
  if (fastq.size() == 2){
    m = 1
  }
"""
danpos.py dpeak -m ${m} \
  ${params.dpeak} \
  -b ${wig_ip_a}:${wig_wce_a},${wig_ip_b}:${wig_wce_b} \
  -o ${file_prefix} \
  ${wig_ip_a}:${wig_ip_b}
"""
}