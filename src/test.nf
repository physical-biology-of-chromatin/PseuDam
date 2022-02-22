
nextflow.enable.dsl=2

include { gatc_finder } from "./nf_modules/gatc_finder/main.nf"
include { htseq_count } from "./nf_modules/htseq_count/main.nf"

workflow {
    gatc_finder(params.fasta)
    

}
