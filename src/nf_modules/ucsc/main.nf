version = "407"
container_url = "lbmc/ucsc:${version}"

params.bedgraph_to_bigwig = ""
params.bedgraph_to_bigwig_out = ""
process bedgraph_to_bigwig {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${file_id}"
  if (params.bedgraph_to_bigwig_out != "") {
    publishDir "results/${params.bedgraph_to_bigwig_out}", mode: 'copy'
  }

  input:
  tuple val(file_id) path(bg)
  tuple val(file_id) path(bed)

  output:
  tuple val(file_id), path("*.bw"), emit: bw

  script:
"""
LC_COLLATE=C
# transform bed file of start-stop chromosome size to stop chromosome size
awk -v OFS="\\t" '{print \$1, \$3}' ${bed} > chromsize.txt

sort -T ./ -k1,1 -k2,2n ${bg} > \
  bedGraphToBigWig ${params.bedgraph_to_bigwig} - \
    chromsize.txt \
    ${bg.simpleName}_norm.bw
"""
}


params.bigwig_to_wig = ""
params.bigwig_to_wig_out = ""
process bigwig_to_wig {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${file_id}"
  if (params.bigwig_to_wig_out != "") {
    publishDir "results/${params.bigwig_to_wig_out}", mode: 'copy'
  }

  input:
  tuple val(file_id) path(bw)

  output:
  tuple val(file_id), path("*.wig"), emit: wig

  script:
"""
bigWigToBedGraph ${bw} ${bw.simpleName}.bg
bedgraph_to_wig.pl --bedgraph ${bw.simpleName}.bg --wig ${bw.simpleName}.wig --step 10
"""
}

params.bigwig2_to_wig2 = ""
params.bigwig2_to_wig2_out = ""
process bigwig2_to_wig2 {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${file_id}"
  if (params.bigwig_to_wig_out != "") {
    publishDir "results/${params.bigwig_to_wig_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(bw_a), path(bw_b)

  output:
  tuple val(file_id), path("${bw_a.simpleName}.wig"), path("${bw_b.simpleName}.wig"), emit: wig

  script:
"""
bigWigToBedGraph ${bw_a} ${bw_a.simpleName}.bg
bedgraph_to_wig.pl --bedgraph ${bw_a.simpleName}.bg --wig ${bw_a.simpleName}.wig --step 10
bigWigToBedGraph ${bw_b} ${bw_b.simpleName}.bg
bedgraph_to_wig.pl --bedgraph ${bw_b.simpleName}.bg --wig ${bw_b.simpleName}.wig --step 10
"""
}

params.wig_to_bigwig = ""
params.wig_to_bigwig_out = ""
process wig_to_bigwig {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${file_id}"
  if (params.bigwig_to_wig_out != "") {
    publishDir "results/${params.bigwig_to_wig_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(w)

  output:
  tuple val(file_id), path("${bw.simpleName}.bw"), emit: bw

  script:
"""
wigToBigWig ${w} ${w.simpleName}.bw
"""
}

params.wig2_to_bigwig2 = ""
params.wig2_to_bigwig2_out = ""
process wig2_to_bigwig2 {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${file_id}"
  if (params.bigwig_to_wig_out != "") {
    publishDir "results/${params.bigwig_to_wig_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(w_a), path(w_b)

  output:
  tuple val(file_id), path("${bw_a.simpleName}.bw"), path("${bw_b.simpleName}.bw"), emit: bw

  script:
"""
wigToBigWig ${w_a} ${w_a.simpleName}.bw
wigToBigWig ${w_b} ${w_b.simpleName}.bw
"""
}