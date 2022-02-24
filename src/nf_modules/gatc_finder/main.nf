nextflow.enable.dsl=2
container_url =  "nathanlecouvreur/gatc_finder"


params.gatc_finder_out = ""

process gatc_finder {
    container = "${container_url}"

  if (params.gatc_finder_out != "") {
    publishDir "results/${params.gatc_finder_out}", mode: 'copy'
  }

    input:
        tuple val(file_id), path(genome)

    output:
        tuple val(file_id), path("*.bed"), emit: bed
        tuple val(file_id), path("*.gff"), emit: gff

    script:
"""
gatc_finder.py --genome ${genome}
"""
}
