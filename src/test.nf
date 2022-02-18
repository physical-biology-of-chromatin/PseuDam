
nextflow.enable.dsl=2

include { gatc_finder } from "./nf_modules/gatc_finder/main.nf"

workflow {
    gatc_finder(params.fasta)
    gatc_finder.out.gatc_sites.view
}
