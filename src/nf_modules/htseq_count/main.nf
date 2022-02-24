nextflow.enable.dsl=2
container_url = "lbmc/htseq:0.13.5"

params.htseq_count_out = ""

process htseq_count {
    container = "${container_url}"

  if (params.htseq_count_out != "") {
    publishDir "results/${params.htseq_count_out}", mode: 'copy'
  }

    input:

        tuple val(file_id), path(bam), path(index)
        tuple val(file_id), path(sites)

    output:

        tuple val(file_id), path("*.count"), emit: count

    script:

"""
htseq-count -f bam -r pos ${bam} ${sites} \
 > results.count
"""
}