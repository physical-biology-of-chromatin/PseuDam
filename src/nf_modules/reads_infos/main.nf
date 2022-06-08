nextflow.enable.dsl=2
container_url =  "nathanlecouvreur/reads_infos"


params.reads_infos_out = ""
params.reads_infos = ""

process reads_infos {
    container = "${container_url}"

  if (params.reads_infos_out != "") {
    publishDir "results/${params.reads_infos_out}", mode: 'copy'
  }

    input:
        tuple val(file_id), path(fastq_file)

    output:
        env(mean), emit: mean
        env(std), emit: std

    script:
"""
test="hello"
out='reads_infos.py --fastq ${fastq_file} ${params.reads_infos}'
echo $test
mean=$out[0]
std=$out[1]
"""
}