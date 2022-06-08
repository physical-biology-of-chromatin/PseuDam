nextflow.enable.dsl=2

params.single_end = ""
params.fasta = "/datas/nathan/vscode_nextflow/nextflow-nathan/data/genome/S288C_reference_sequence_R64-3-1_20210421.fsa"
params.fastq = "/datas/nathan/vscode_nextflow/nextflow-nathan/data/reads/Dam_ID/split_DamC/*{1,2}.fq"
params.bed = ""


include { fastp                   } from "./nf_modules/fastp/main.nf"       
include { gatc_finder             } from "./nf_modules/gatc_finder/main.nf" addParams(gatc_finder_out: "Dam_ID/DamC_split/sites_test/", gatc_finder: "--overlap")
include { index_fasta; mapping_fastq} from "./nf_modules/kallisto/main.nf"    addParams(mapping_fastq_out: "Dam_ID/DamC_split/counts/", mapping_fastq: "--bias --bootstrap-samples 100 ${params.single_end}")
include { multiqc                 } from "./nf_modules/multiqc/main.nf"     addParams(multiqc_out: "Dam_ID/DamC_split/reports/")
include { fasta_from_bed          } from "./nf_modules/bedtools/main.nf"    addParams(fasta_from_bed_out: "Dam_ID/DamC_split/fasta/", fasta_from_bed: "-name")




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

    fasta_from_bed(fasta_files,
                   gatc_finder.out.bed)

    index_fasta(fasta_from_bed.out.fasta)

    fastp(fastq_files)

    mapping_fastq(index_fasta.out.index.collect(),
                  fastp.out.fastq)

    report_mapping = mapping_fastq.out.report

    report_fastp = fastp.out.report

    multiqc(report_mapping.mix(report_fastp))

}