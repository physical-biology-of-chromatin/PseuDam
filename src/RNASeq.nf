nextflow.enable.dsl=2


include { fastp } from "./nf_modules/fastp/main.nf"

params.fastq = "data/fastq/*_{1,2}.fastq"
channel
  .fromFilePairs( params.fastq, size: -1)
  .set { fastq_files }


workflow {
    fastp(fastq_files)
}
