Channel
  .fromPath( "data/tiny_dataset/fasta/*.fasta" )
  .set { fasta_file }

process sample_fasta {
  publishDir "results/sampling/", mode: 'copy'

  input:
file fasta from fasta_file

  output:
file "sample.fasta" into fasta_sample

  script:
"""
head ${fasta} > sample.fasta
"""
}
