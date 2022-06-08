nextflow.enable.dsl=2
container_url =  "nathanlecouvreur/fastq_splitter"


params.fastq_splitter_out = ""
params.fastq_splitter = ""


workflow split {
  take:
  reads

  main:
    split_fastq(reads)
    group_fastq(split_fastq.out.fastq_folder)
    group_fastq.out.fastq
      .map{it -> it[1]}
      .flatten()
      .collate(2)
      .map{ it -> [it[0].simpleName - ~/_{0,1}R[12]/, it]}
      .set{ splited_fastq }

  emit:
    fatq = splited_fastq
}


process split_fastq {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  if (params.fastq_splitter_out != "") {
    publishDir "results/${params.fastq_splitter_out}", mode: 'copy'
  }

    input:
      tuple val(file_id), path(fastq)


    output:
      tuple val(file_id), path("split"), emit: fastq_folder

    script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }

  if (reads.size() == 2)
  """
  fastq_splitter ${params.fastq_splitter} -n 2 ${reads[0]},${reads[1]}
  mv *.fq split

  """


process group_fastq {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  if (params.fastq_splitter_out != "") {
    publishDir "results/${params.split_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(reads_folder)

  output:
  tuple val(file_id), path("results/*"), emit: fastq

  script:
"""
mkdir -p results/
find split/ -type "f" | \
  grep -v "unassigned" | \
  sed -E "s|(split/(.*)/(.*))|\\1 \\2_\\3|g" |
  awk '{system("mv "\$1" results/"\$2)}'
"""
}




}