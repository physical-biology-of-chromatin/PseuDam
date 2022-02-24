nextflow.enable.dsl=2
/*
./nextflow src/Dam_ID_analysis.nf -profile docker --fasta data/genome/dm6.fasta --fastq data/reads/data_R.fastq
*/


include { fastp } from "./nf_modules/fastp/main.nf" 
include { index_fasta ; mapping_fastq } from "./nf_modules/bowtie2/main.nf" 
include { index_bam ; sort_bam} from "./nf_modules/samtools/main.nf" 
include { gatc_finder } from "./nf_modules/gatc_finder/main.nf"
include { multiqc } from "./nf_modules/multiqc/main.nf" addParams(multiqc_out: "mapping/")
include { htseq_count } from "./nf_modules/htseq_count/main.nf"


params.fasta = "data/genome/*_G.fasta"
params.fastq = "data/reads/*_R{1,2}.fastq"


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

    gatc_finder(fasta_files)

    mapping_fastq(index_fasta.out.index.collect(), 
                  fastp.out.fastq)

    sort_bam(mapping_fastq.out.bam)

    index_bam(sort_bam.out.bam)

    report_mapping = mapping_fastq.out.report

    report_fastp = fastp.out.report

    multiqc(report_mapping.mix(report_fastp))

    htseq_count(index_bam.out.bam_idx, 
                gatc_finder.out.gff) 

    
}

