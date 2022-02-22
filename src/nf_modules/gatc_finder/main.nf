nextflow.enable.dsl=2
container_url =  "nathanlecouvreur/gatc_finder"

process gatc_finder {
    container = "${container_url}"


    input:
        file(genome)

    output:
        path "*.bed", emit: bed
        path "*.gff", emit: gff

    script:

"""
gatc_finder.py --genome ${genome}
"""
}
