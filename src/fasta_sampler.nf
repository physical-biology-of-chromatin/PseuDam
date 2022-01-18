nextflow.enable.dsl=2

channel
  .fromPath( " data/tiny_dataset/fasta/*.fasta" )
  .set { fasta_file }

workflow{
  sample_fasta(fasta_file)
}
