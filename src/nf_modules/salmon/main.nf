nextflow.enable.dsl=2
version = "0.8.2"
container_url = "combinelab/salmon:latest"

params.index_fasta = ""
params.index_fasta_out = ""
process index_fasta {
  container = "${container_url}"
  label "big_mem_mono_cpus"

  if (params.index_fasta_out != "") {
    publishDir "results/${params.index_fasta_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(fasta)

  output:
  tuple val(file_id), path("*.index"), emit: index

  script:

  """
  salmon index -t ${fasta} -i ${fasta.baseName}.index ${params.index_fasta}
  """
}


params.mapping_fastq = ""
params.mapping_fastq_out = ""
process mapping_fastq {
    container = "${container_url}"
    label = "big_mem_mono_cpus"

    if (params.mapping_fastq_out != "") {
        publishDir "results/${params.mapping_fastq_out}", mode: "copy"
    }

    input:
        tuple val(index_id), path(index)
        tuple val(file_id), path(reads)

    output:
        tuple val(file_id), path("${file_prefix}"), emit: counts
        tuple val(file_id), path("*_report.txt"), emit: report

    script:
    if (file_id instanceof List){
        file_prefix = file_id[0]
    } 
    else {
        file_prefix = file_id
    }

    if (reads.size() == 2)
    """
    mkdir ${file_prefix}
    salmon quant -i ${index} \
    ${params.mapping_fastq} -o ${file_prefix} \
    -1 ${reads[0]} -2 ${reads[1]} &> ${file_prefix}_salmon_mapping_report.txt
    """
    else
    """
    mkdir ${file_prefix}
    salmon quant -i ${index} \
    ${params.mapping_fastq} -o ${file_prefix} \
    -r ${reads[0]} &> ${file_prefix}_salmon_mapping_report.txt
    """
}