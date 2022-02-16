
nextflow.enable.dsl=2

include { gatc_finder } from "./nextflow-nathan/src/nf_modules/gatc_finder/main.nf"


params.fasta = "data/genome/*_G.fasta"

workflow {
    gatc_finder(params.fasta)
    gatc_finder.out.gatc_sites.view

}
