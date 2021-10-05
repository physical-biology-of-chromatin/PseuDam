version = "v2.2.2_cv3"
container_url = "biocontainers/danpos:${version}"

include {
  bigwig_to_wig
} from "./../ucsc/main.nf"

params.dpos = "--smooth_width 0 -n N "
params.dpos_out = ""

process dpos_bam {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.dpos_out != "") {
    publishDir "results/${params.dpos_out}", mode: 'copy', pattern: "${file_prefix}/*"
  }

  input:
    tuple val(file_id), path(fastq)
    tuple val(file_id), path(bam_ip)
    tuple val(file_id), path(bam_wce)

  output:
    tuple file_id, "${file_prefix}/pooled/*.wig" into wig
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
danpos.py dpos -m ${m}
  ${params.dpos} \
  -b ${bam_wce} \
  -o ${file_prefix} \
  ${bam_ip}
"""
}

workflow dpos_bw {
  take:
    fastq
    bw_ip
    bw_wce
  main:
    dpos_wig(fastq, bigwig_to_wig(bw_ip), bigwig_to_wig(bw_wce))
  emit:
  wig = dpos_wig.out.wig
  folder = dpos_wig.out.folder
}

process dpos_wig {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.dpos_out != "") {
    publishDir "results/${params.dpos_out}", mode: 'copy', pattern: "${file_prefix}/*"
  }

  input:
    tuple val(file_id), path(fastq)
    tuple val(file_id), path(wig_ip)
    tuple val(file_id), path(wig_wce)

  output:
  tuple file_id, "${file_prefix}/pooled/*.wig" into wig
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
danpos.py dpos -m ${m}
  ${params.dpos} \
  -b ${wig_wce} \
  -o ${file_prefix} \
  ${wig_ip}
"""
}

process dpos_wigvswig {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.dpos_out != "") {
    publishDir "results/${params.dpos_out}", mode: 'copy', pattern: "${file_prefix}/*"
  }

  input:
    tuple val(file_id), path(fastq)
    tuple val(file_id), path(wig_ip_a)
    tuple val(file_id), path(wig_wce_a)
    tuple val(file_id), path(wig_ip_b)
    tuple val(file_id), path(wig_wce_b)

  output:
  tuple file_id, "${file_prefix}/pooled/*.wig" into wig
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
danpos.py dpos -m ${m}
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
    publishDir "results/${params.dpeak_out}", mode: 'copy', pattern: "${file_prefix}/*"
  }

  input:
    tuple val(file_id), path(fastq)
    tuple val(file_id), path(bam_ip)
    tuple val(file_id), path(bam_wce)

  output:
    tuple file_id, "${file_prefix}/pooled/*.wig" into wig
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
danpos.py dpeak -m ${m}
  ${params.dpeak} \
  -b ${bam_wce} \
  -o ${file_prefix} \
  ${bam_ip}
"""
}

workflow dpeak_bw {
  take:
    fastq
    bw_ip
    bw_wce
  main:
    dpeak_wig(fastq, bigwig_to_wig(bw_ip), bigwig_to_wig(bw_wce))
  emit:
  wig = dpeak_wig.out.wig
  folder = dpeak_wig.out.folder
}


process dpeak_wig {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.dpeak_out != "") {
    publishDir "results/${params.dpeak_out}", mode: 'copy', pattern: "${file_prefix}/*"
  }

  input:
    tuple val(file_id), path(fastq)
    tuple val(file_id), path(wig_ip)
    tuple val(file_id), path(wig_wce)

  output:
  tuple file_id, "${file_prefix}/pooled/*.wig" into wig
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
danpos.py dpeak -m ${m}
  ${params.dpeak} \
  -b ${wig_wce} \
  -o ${file_prefix} \
  ${wig_ip}
"""
}

workflow dpeak_bwvsbw {
  take:
    fastq
    bw_ip_a
    bw_wce_a
    bw_ip_b
    bw_wce_b
  main:
    dpeak_wigvswig(
      fastq,
      bigwig_to_wig(bw_ip_a),
      bigwig_to_wig(bw_wce_a),
      bigwig_to_wig(bw_ip_b),
      bigwig_to_wig(bw_wce_b)
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
    publishDir "results/${params.dpeak_out}", mode: 'copy', pattern: "${file_prefix}/*"
  }

  input:
    tuple val(file_id), path(fastq)
    tuple val(file_id), path(wig_ip_a)
    tuple val(file_id), path(wig_wce_a)
    tuple val(file_id), path(wig_ip_b)
    tuple val(file_id), path(wig_wce_b)

  output:
  tuple file_id, "${file_prefix}/pooled/*.wig" into wig
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
danpos.py dpeak -m ${m}
  ${params.dpeak} \
  -b ${wig_ip_a}:${wig_wce_a},${wig_ip_b}:${wig_wce_b} \
  -o ${file_prefix} \
  ${wig_ip_a}:${wig_ip_b}
"""
}