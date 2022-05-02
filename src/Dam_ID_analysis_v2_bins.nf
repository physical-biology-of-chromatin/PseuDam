nextflow.enable.dsl=2

params.bin_dir = ""
params.mean_length = ""
params.std_length = ""
params.fasta = "data/genome/GCA_000002985.3.fasta"
params.fastq = "/datas/nathan/vscode_nextflow/nextflow-nathan/data/reads/Dam_ID/c_elegans/dpy-7/*_L2_*.fastq"
params.bed = ""


include { fastp                   } from "./nf_modules/fastp/main.nf"       
include { gatc_finder             } from "./nf_modules/gatc_finder/main.nf" addParams(gatc_finder_out: "Dam_ID/c_elegans/binning_test/${params.bin_dir}/sites_test/", gatc_finder: "--overlap")
include { index_fasta; mapping_fastq} from "./nf_modules/kallisto/main.nf"    addParams(mapping_fastq_out: "Dam_ID/c_elegans/binning_test/${params.bin_dir}/counts/", mapping_fastq: "--bias --bootstrap-samples 100 -l ${params.mean_length} -s ${params.std_length}")
include { multiqc                 } from "./nf_modules/multiqc/main.nf"     addParams(multiqc_out: "Dam_ID/c_elegans/binning_test/${params.bin_dir}/reports/")
include { fasta_from_bed          } from "./nf_modules/bedtools/main.nf"    addParams(fasta_from_bed_out: "Dam_ID/c_elegans/binning_test/${params.bin_dir}/fasta/", fasta_from_bed: "-name")


channel
    .fromPath(params.bed)
    .ifEmpty { error "Cannot find any bed file matching: ${params.fasta}" }
    .map { it -> [it.simpleName, it]}
    .set{bed}

channel
    .fromPath(params.fasta)
    .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
    .map { it -> [it.simpleName, it]}
    .set {fasta_files}

channel
    .fromFilePairs(params.fastq, size: -1)
    .set {fastq_files}


workflow {

    fasta_from_bed(fasta_files,
                   bed)

    index_fasta(fasta_from_bed.out.fasta)

    fastp(fastq_files)

    mapping_fastq(index_fasta.out.index.collect(),
                  fastp.out.fastq)

    report_mapping = mapping_fastq.out.report

    report_fastp = fastp.out.report

    multiqc(report_mapping.mix(report_fastp))

}