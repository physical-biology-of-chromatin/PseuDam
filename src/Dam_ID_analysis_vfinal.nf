nextflow.enable.dsl=2


def help_message() {

    log.info"""
    =========================================
                 Dam_ID_analysis
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    ./nextflow main.nf -profile docker --fastq '/data/*_{R1,R2}*.fastq' --outdir '/project/'

    Required arguments:
         -profile                      Configuration profile to use. <base, docker>
         --fastq                       Directory pattern for fastq files: /project/*{R1,R2}*.fastq
         --fasta                       Fasta file of the reference genome to use

    Single end options:
         --single_end                  Sets the pipeline to accept single end datas. Requires to input --mean_length and --std_length
         --mean_length                 Mean fragment length of the samples (if not accessible, use the reads' mean length)
         --std_length                  Std of the fragment length of the samples (if not accessible, use the reads' length std or 10 if it is 0)

    Binning option:
        --bed                          Provides a .bed file with sites instead of processing it from the reference genome with GATC_finder   
        
    QC Options:
        --skipFastp                    Skip running fastp
    """.stripIndent()
}


params.help = false

if (params.help != false){
    help_message()
    exit 0
}


params.outdir = ""

params.fasta = "/datas/nathan/vscode_nextflow/nextflow-nathan/data/genome/S288C_reference_sequence_R64-3-1_20210421_12less.fsa"
params.fastq = "/datas/nathan/vscode_nextflow/nextflow-nathan/data/reads/Dam_ID/split_DamC/*{1,2}.fq"

params.skipFastp = false



params.bed = ""

if (params.bed != ""){
    channel
        .fromPath(params.bed)
        .ifEmpty { error "Cannot find any bed file matching: ${params.fasta}" }
        .map { it -> [it.simpleName, it]}
        .set{bed}
}

channel
    .fromPath(params.fasta)
    .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
    .map { it -> [it.simpleName, it]}
    .set {fasta_files}

channel
    .fromFilePairs(params.fastq, size: -1)
    .set {fastq_files}



include { fastp                      } from "./nf_modules/fastp/main.nf"       
include { gatc_finder                } from "./nf_modules/gatc_finder/main.nf" addParams(gatc_finder_out: "${params.outdir}/sites_test/", gatc_finder: "--overlap")
include { index_fasta; mapping_fastq } from "./nf_modules/kallisto/main.nf"    addParams(mapping_fastq_out: "${params.outdir}/counts/")
include { multiqc                    } from "./nf_modules/multiqc/main.nf"     addParams(multiqc_out: "${params.outdir}/reports/")
include { fasta_from_bed             } from "./nf_modules/bedtools/main.nf"    addParams(fasta_from_bed_out: "${params.outdir}/fasta/", fasta_from_bed: "-name")



params.single_end = false

params.mean_length = ""
params.std_length = ""


if (params.single_end != false){

    if (params.mean_length == ""){
        error("In case of single end datas, please use --mean_length and --std_length to provide the mean and std of the fragments (if not known use the reads mean length and std)")
    } 
    if (params.std_length == "") {
         error("In case of single end datas, please use --mean_length and --std_length to provide the mean and std of the fragments (if not known use the reads mean length and std)")
    }

    addParams(mapping_fastq: "-l ${params.mean_length} -s ${params.std_length} --bias --bootstrap-samples 100")
}



workflow {


    if (params.bed != ""){

        fasta_from_bed(fasta_files,
                       bed)

    } else {

        gatc_finder(fasta_files)

        fasta_from_bed(fasta, 
                       gatc_finder.out.bed)
    }


    index_fasta(fasta_from_bed.out.fasta)


    if (params.skipFastp != false){

        fastp(fastq_files)

        mapping_fastq(index_fasta.out.index.collect(),
                      fastq_files)

        multiqc(mapping_fastq.out.report)

    } else {
        
        mapping_fastq(index_fasta.out.index.collect(),
                      fastp.out.fastq)

        report_mapping = mapping_fastq.out.report

        report_fastp = fastp.out.report

        multiqc(report_mapping.mix(report_fastp))
    }

}