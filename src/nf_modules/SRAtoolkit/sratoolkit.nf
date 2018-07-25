/*
* sra-tools :

*/

/*                      fastq-dump
* Imputs : srr list
* Outputs : fastq files
*/

params.list_srr = "$baseDir/data/SRR/*.txt"

log.info "downloading list srr : ${params.list_srr}"

Channel
  .fromPath( params.list_srr )
  .ifEmpty { error "Cannot find any bam files matching: ${params.list_srr}" }
  .splitCsv(header: true)
  .set { SRR }

//run is the column name containing SRR ids

process fastq_dump {
  tag {"${x.run}"}
  publishDir "results/download/fastq/${x.run}/", mode: 'copy'

  input:
    val x  from SRR

  output:
    file("*") into fastq

  script:
"""
fastq-dump --split-files --defline-seq '@\$ac_\$si/\$ri' --defline-qual "+"  ${x.run}
if [ -f ${x.run}_1.fastq ]
then
  true
else
  touch ${x.run}.fastq
fi
"""
}
