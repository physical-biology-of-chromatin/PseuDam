nextflow.enable.dsl=2
container_url =  "nathanlecouvreur/gatc_finder"

process gatc_finder {
    container = "${container_url}"


    input:
        path(genome)

    output:
        path "*.bed", emit: gatc_sites

    script:

"""
python3 gatc_finder.py --genome ${genome}
"""
}
