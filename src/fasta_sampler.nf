nextflow.enable.dsl=2

channel
  .fromPath( " data/tiny_dataset/fasta/*.fasta" )
  .set { fasta_file }

process sample_fasta {

  input:
  path fasta
  

  output:
  path "sample.fasta", emit: fasta_sample
  publishDir "results/sampling/", mode: 'copy'
  
  script:
"""
head ${fasta} > ${fasta.simpleName}_sample.fasta
"""
}

workflow{
  sample_fasta(fasta_file)
}
