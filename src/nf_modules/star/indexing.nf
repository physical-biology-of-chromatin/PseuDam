params.fasta = "$baseDir/data/bam/*.fasta"
params.annotation = "$baseDir/data/bam/*.gtf"

log.info "fasta files : ${params.fasta}"

Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
  .set { fasta_file }
Channel
  .fromPath( params.annotation )
  .ifEmpty { error "Cannot find any annotation files matching: ${params.annotation}" }
  .set { annotation_file }

process index_fasta {
  tag "$fasta.baseName"
  publishDir "results/mapping/index/", mode: 'copy'

  input:
    file fasta from fasta_file
    file annotation from annotation_file

  output:
    file "*" into index_files

  script:
"""
STAR --runThreadN ${task.cpus} --runMode genomeGenerate \
--genomeDir ./ \
--genomeFastaFiles ${fasta} \
--sjdbGTFfile ${annotation} \
--genomeSAindexNbases 3 # min(14, log2(GenomeLength)/2 - 1)
"""
}


