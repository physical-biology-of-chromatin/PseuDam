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
  tuple val(file_id) path(bg)

  output:
  tuple val(file_id), path("*.wig"), emit: wig

  script:
"""
bigWigToBedGraph ${params.bigwig_to_wig} \
  ${bg} \
  ${bg.simpleName}.bg
awk '{if(NR>1) {if(\$1!=lastChrom){printf("variableStep chrom=%s\\n",\$1);lastChrom=\$1;}print \$2,\$4}}' ${bg.simpleName}.bg > ${bg.simpleName}.wig
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
bigWigToWig ${params.bigwig_to_wig} \
  ${bw_a} \
  ${bw_a.simpleName}.wig
bigWigToWig ${params.bigwig_to_wig} \
  ${bw_b} \
  ${bw_b.simpleName}.wig
"""
}