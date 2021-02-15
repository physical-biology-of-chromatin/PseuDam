version = "407"
container_url = "lbmc/ucsc:${version}"

process bedgraph_to_bigwig {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${file_id}"

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
  bedGraphToBigWig - \
    chromsize.txt \
    ${bg.simpleName}_norm.bw
"""
}
