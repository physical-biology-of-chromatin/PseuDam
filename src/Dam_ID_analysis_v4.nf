nextflow.enable.dsl=2

include { fastp                             } from "./nf_modules/fastp/main.nf"
include { gatc_finder                       } from "./nf_modules/gatc_finder/main.nf" addParams(gatc_finder: "--salmon", gatc_finder_out: "sites/salmon/")
include { index_fasta; mapping_fastq        } from "./nf_modules/salmon/main.nf"      addParams(mapping_fastq_out: "salmon/")
include { multiqc                           } from "./nf_modules/multiqc/main.nf"     addParams(multiqc_out: "mapping/salmon/")
include { fasta_from_bed                    } from "./nf_modules/bedtools/main.nf"


params.fasta = "data/genome/S288C_reference_sequence_R64-3-1_20210421_12less.fsa"
params.fastq = "data/reads/Dam_ID/N_3RG1D2_{1,2}.fq"


channel
    .fromPath(params.fasta)
    .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
    .map { it -> [it.simpleName, it]}
    .set {fasta_files}

channel
    .fromFilePairs(params.fastq, size: -1)
    .set {fastq_files}


workflow {

    gatc_finder(fasta_files)

    fastp(fastq_files)

    fasta_from_bed(fasta_files,
                   gatc_finder.out.bed)

    index_fasta(fasta_from_bed.out.fasta)

    mapping_fastq(index_fasta.out.index.collect(),
                  fastp.out.fastq)

    report_mapping = mapping_fastq.out.report

    report_fastp = fastp.out.report

    multiqc(report_mapping.mix(report_fastp))    
}