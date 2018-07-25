params.fastq = "$baseDir/data/fastq/*.fastq"

log.info "fastq files : ${params.fastq}"

Channel
  .fromPath( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

process fastqc_fastq {
  tag "$reads.baseName"
  publishDir "results/fastq/fastqc/", mode: 'copy'
  cpus = 1

  input:
    file reads from fastq_files

  output:
    file "*.{zip,html}" into fastqc_repport

  script:
"""
fastqc --quiet --threads ${task.cpus} --format fastq --outdir ./ ${reads}
"""
}

process multiqc {
  tag "$repport[0].baseName"
  publishDir "results/fastq/multiqc/", mode: 'copy'
  cpus = 1

  input:
    file repport from fastqc_repport.collect()

  output:
    file "*multiqc_*" into multiqc_report

  script:
"""
multiqc -f .
"""
}

