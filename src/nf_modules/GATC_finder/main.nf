params.genome = ""


channel
    .fromPath(params.genome)
    .set(genome)

process GATC_finder {
    container = "/home/nathan/projects/vscode_nextflow/nextflow-nathan/src/.docker_modules/GATC_finder"
    label "?"
    tag "?"
}

    input:
        file fasta from genome

    output:
        file sites into gatc_sites

"""
docker run --volume ${fasta}:/genome.fa gatc_finder \

    gatc_finder genome.fa
"""