nextflow.enable.dsl=2
/*
./nextflow src/Dam_ID_analysis.nf -profile docker --reads "data/reads/data.fastq --genome "data/genome/dm6.fa
*/





/*========================= modules import ================================*/

include { fastp } from "./nf_modules/fastp/main.nf"

include { index_fasta; mapping_fastq } from "./nf_modules/bowtie2/main.nf"


params.reads = ""
params.genome = ""

log.info "reads files : ${params.reads}"
log.info "genome file : ${params.genome}"


/*
channel
    .fromPath( params.bed )
    .ifEmpty { error "Cannot find any bed files matching: ${params.bed}" }
    .map { it -> [it.simpleName, it]}
    .set { bed_files }
*/

channel
    .fromPath( params.reads )
    .set {reads}


/*================================ workflow ================================*/

workflow {
    fastp(reads)
    //mapping
    index_fasta(reads)
    mapping_fastq(index_fasta.out.out.index.collect(), 
                  reads)
}

workflow.onComplete {
    println ( workflow.success "Au moins ça a pas planté, mais ça a surement pas marché quand même" )
}
