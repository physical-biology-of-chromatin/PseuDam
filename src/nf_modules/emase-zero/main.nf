version = "0.3.1"
container_url = "lbmc/emase-zero:${version}"

include { tr2g } from "./../kb/main.nf"


params.count = "-m 2"
params.count_out = ""
workflow count {
  take:
    bin
    transcript_length
    gtf

  main:
    tr2g(gtf, channel.of(["NO T2G", ""]))
    emase(bin, transcript_length, tr2g.out.t2g)

  emit:
    count: emase.count
}

process emase {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.count_out != "") {
    publishDir "results/${params.count_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bin)
    tuple val(transcript_length_id), path(transcript_length)
    tuple val(transcript_to_gene_id), path(transcript_to_gene)

  output:
    tuple val(file_id), path("${bin.simpleName}.quantified"), emit: count

  script:
"""
emase-zero ${params.count} \
  -b ${bin} \
  -o ${bin.simpleName}.quantified \
  -l ${transcript_length} \
  -g ${transcript_to_gene}
"""
}