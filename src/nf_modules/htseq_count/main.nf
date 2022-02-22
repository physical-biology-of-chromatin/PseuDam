nextflow.enable.dsl=2
container_url = "dmccloskey/htseq-count"

process htseq_count {
    container = "${container_url}"


    input:

        path(reads)
        path(sites)

    output:

        path "*.csv", emit: count

    script:

"""
htseq-count -f bam ${reads} ${sites}
"""
}