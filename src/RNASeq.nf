nextflow.enable.dsl=2

include { fastp } from "./nf_modules/fastp/main.nf"

channel
  .fromFilePairs( "data/tiny_dataset/fastq/*_R{1,2}.fastq", size: -1)
  .set { fastq_files }

workflow {
    fastp(fastp_files)
}
