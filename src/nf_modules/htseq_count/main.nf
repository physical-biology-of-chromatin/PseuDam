nextflow.enable.dsl=2
container_url = "lbmc/htseq:0.13.5"

params.htseq_count = ""
params.htseq_count_out = ""

process htseq_count {
    container = "${container_url}"

  if (params.htseq_count_out != "") {
    publishDir "results/${params.htseq_count_out}", mode: 'copy'
  }

    input:

        tuple val(bam_id), path(bam), path(index)
        tuple val(sites_id), path(sites)

    output:

        tuple val(bam_id), path("*"), emit: count

    script:

"""
htseq-count -f bam -r pos ${bam} ${sites} ${params.htseq_count} \
> ${bam_id}.txt
"""
}