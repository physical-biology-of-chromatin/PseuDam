nextflow.enable.dsl=2
/*
./nextflow src/Dam_ID_analysis.nf -profile docker --fasta <reference genome> --fastq <genome>
*/

include { fastp                             } from "./nf_modules/fastp/main.nf" 
include { index_fasta ; mapping_fastq       } from "./nf_modules/bowtie2/main.nf" 
include { index_bam ; sort_bam              } from "./nf_modules/samtools/main.nf"    
include { gatc_finder                       } from "./nf_modules/gatc_finder/main.nf" addParams(gatc_finder_out: "sites/")
include { multiqc                           } from "./nf_modules/multiqc/main.nf"     addParams(multiqc_out: "mapping/")
include { bed_to_gff                        } from "./nf_modules/gffread/main.nf"     
include { htseq_count                       } from "./nf_modules/htseq_count/main.nf" addParams(htseq_count_out: "htseq_count/")


/* Input parameters */

params.fasta = "data/genome/S288C_reference_sequence_R64-3-1_20210421.fsa"
params.fastq = "data/reads/Dam_ID/N_3RG3D2_EKDL220001512-1a_H2C2YDSX3_L3_{1,2}.fq"




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

    /* Localisation of all the gatc sites in the reference genome and creation of the bins */
    gatc_finder(fasta_files)

    /* Quality control and adaptator trimming */
    fastp(fastq_files)

    /* Indexing of the reference genome */
    index_fasta(fasta_files)

    /* mapping of the reads to the indexed genome */
    mapping_fastq(index_fasta.out.index.collect(),
                  fastp.out.fastq)

    /* sorting of the mapped reads */
    sort_bam(mapping_fastq.out.bam)


    /* indexing of the mapped reads */
    index_bam(sort_bam.out.bam)


    /* creation of a global report of the qc and mapping */
    report_mapping = mapping_fastq.out.report


    report_fastp = fastp.out.report


    multiqc(report_mapping.mix(report_fastp))

    htseq_count(index_bam.out.bam_idx,
                gatc_finder.out.gff)
    
}

