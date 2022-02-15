container_url =  "nathanlecouvreur/gatc_finder"
params.fasta = ""


channel
    .fromPath(params.fasta)
    .view()
    .set{genome}


process GATC_finder {
    container = "${container_url}"


    input:
        file genome from genome

    output:
        file sites into gatc_sites


"""
docker run --volume ${genome}:/genome.fa gatc_finder \
    gatc_finder genome.fa
"""
}
