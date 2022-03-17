nextflow.enable.dsl=2
/*
./nextflow src/Dam_ID_analysis.nf -profile docker --fasta <reference genome> --fastq <genome>
*/

include { fastp                             } from "./nf_modules/fastp/main.nf" 
include { index_fasta ; mapping_fastq       } from "./nf_modules/bowtie2/main.nf" 
include { index_bam ; sort_bam              } from "./nf_modules/samtools/main.nf"    
include { gatc_finder                       } from "./nf_modules/gatc_finder/main.nf"
include { multiqc                           } from "./nf_modules/multiqc/main.nf"     addParams(multiqc_out: "mapping/")
include { bam_to_bed ; intersect            } from "./nf_modules/bedtools/main.nf"    addParams(bam_to_bed_out: "raw_bed/", intersect_out: "intersect/", bam_to_bed: "-bedpe")
include { coverage as coverage_base         } from "./nf_modules/bedtools/main.nf"    addParams(coverage_out: "coverage_base/", coverage: "-d")
include { coverage as coverage_frag         } from "./nf_modules/bedtools/main.nf"    addParams(coverage_out: "coverage_frag/")
include { bed_to_gff                        } from "./nf_modules/gffread/main.nf"     



/* Input parameters */

params.fasta = "data/genome/S288C_reference_sequence_R64-3-1_20210421.fsa"
params.fastq = "data/reads/Dam_ID/*_{1,2}.fq"




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

    /* Quality control and adaptator trimming */
    fastp(fastq_files)


    /* Localisation of all the gatc sites in the reference genome and creation of the bins */
    gatc_finder(fasta_files)



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


    /* conversion of the bam files to bed files
    bam_to_bed(index_bam.out.bam_idx)
    */

    /* selection of the reads entirely in a single bin */
    intersect(index_bam.out.bam_idx, 
              gatc_finder.out.bed.collect())


    /* calculation of the coverage of each fragment */
    coverage_frag(intersect.out.intersect,
                  gatc_finder.out.bed.collect())

    /* calculation of the coverage for each base in each fragment */
    coverage_base(intersect.out.intersect,
                  gatc_finder.out.bed.collect())
}

