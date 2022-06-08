
nextflow.enable.dsl=2

include {gatc_finder; frag_length} from "/datas/nathan/vscode_nextflow/nextflow-nathan/src/nf_modules/gatc_finder/main.nf" 

params.fasta = "/datas/nathan/vscode_nextflow/nextflow-nathan/data/genome/GCA_000002985.3.fasta"

channel
    .fromFilePairs(params.fasta, size: -1)
    .set {fasta_files}


workflow {
    gatc_finder(fasta_files)
    frag_length(gatc_finder.out.bed)
    frag_length.out.mean.view()
    frag_length.out.std.view()
}
