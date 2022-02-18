nextflow.enable.dsl=2
/*
./nextflow src/Dam_ID_analysis.nf -profile docker --reads "data/reads/data.fastq --genome "data/genome/dm6.fa
*/


include { fastp } from "./nf_modules/fastp/main.nf" 
include { index_fasta ; mapping_fastq } from "./nf_modules/bowtie2/main.nf" 
include { index_bam ; sort_bam} from "./nf_modules/samtools/main.nf" 
include { gatc_finder } from "./nf_modules/gatc_finder/main.nf"
include { multiqc } from "./nf_modules/multiqc/main.nf" addParams(multiqc_out: "mapping/")

params.fasta = "data/genome/*_G.fasta"
params.fastq = "data/reads/*_R.fastq"


channel
    .fromPath(params.fasta)
    .ifEmpty { error "Cannot find any fasta files matching: ${params.genome}" }
    .map { it -> [it.simpleName, it]}
    .set {fasta_files}

channel
    .fromFilePairs(params.fastq, size: -1)
    .set {fastq_files}
 
/*================================ workflow ================================*/

workflow {
    
    fastp(fastq_files)

    index_fasta(fasta_files)

    mapping_fastq(index_fasta.out.index.collect(), 
                  fastp.out.fastq)

    sort_bam(mapping_fastq.out.bam)

    index_bam(sort_bam.out.bam)

    multiqc(index_bam.out.bam_idx)

}

