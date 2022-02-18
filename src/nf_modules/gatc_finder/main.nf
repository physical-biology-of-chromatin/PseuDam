nextflow.enable.dsl=2
container_url =  "nathanlecouvreur/gatc_finder"

process gatc_finder {
    container = "${container_url}"


    input:
        file(genome)

    output:
        path "*.bed", emit: gatc_sites

    script:
"""
gatc_finder.py --genome ${genome}
"""
}
