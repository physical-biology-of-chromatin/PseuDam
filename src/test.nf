
nextflow.enable.dsl=2

include {bed_to_gff} from "/datas/nathan/vscode_nextflow/nextflow-nathan/src/nf_modules/gffread/main.nf" 

params.bed = "results/GATC/*_new.bed"

channel
    .fromPath(params.bed)
    .set{bed}


workflow {
    bed_to_gff(bed)
    bed_to_gff.out.view()

}
