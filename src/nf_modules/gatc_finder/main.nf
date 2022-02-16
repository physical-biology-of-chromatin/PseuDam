nextflow.enable.dsl=2
container_url =  "nathanlecouvreur/gatc_finder"
params.fasta = ""


channel
    .fromPath(params.fasta)
    .view()
    .set{genome}


process gatc_finder {
    container = "${container_url}"


    input:
        path(genome)

    output:
        path "*.bed", emit: gatc_sites

    script:
"""
docker run --volume ${genome}:/genome.fa gatc_finder \
    gatc_finder genome.fa
    
"""
}
