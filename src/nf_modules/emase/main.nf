version = "0.10.16"
container_url = "lbmc/emase:${version}"

params.diploid_genome = "-x"
params.diploid_genome_out = "-x"
process diploid_genome {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.diploid_genome_out != "") {
    publishDir "results/${params.diploid_genome_out}", mode: 'copy'
  }

  input:
    tuple val(genome_a), path(fasta_a)
    tuple val(genome_b), path(fasta_b)

  output:
    tuple val(file_id), path(".fa"), emit: index

  script:
"""
prepare-emase -G ${fasta_a},${fasta_b} -s ${genome_a},${genome_b} ${params.diploid_genome} 
"""
}