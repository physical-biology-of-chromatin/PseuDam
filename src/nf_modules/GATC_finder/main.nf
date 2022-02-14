container

params.genome = ""
params.out_file = ""

process GATC_finder {
    container = "/home/nathan/projects/vscode_nextflow/nextflow-nathan/src/.docker_modules/GATC_finder"
    label "?"
    tag "?"
}

    input:
        val params.genome
        val params.out_file

    output:
        file "sites.bed"

"""
gatc_finder ${params.genome} ${params.out_file}
"""