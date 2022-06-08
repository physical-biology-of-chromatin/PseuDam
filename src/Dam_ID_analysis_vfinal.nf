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
         --fastq                       Directory pattern for fastq files: 
                                       /project/*{R1,R2}*.fastq
         --fasta                       Fasta file of the reference genome to use
    
    Additionnal arguments :
         --outdir                      Specifiest the directory to save the files in within /results

    Single end options:
         --single_end                  Sets the pipeline to accept single end datas.
                                       Requires to input --mean_length and --std_length
         --mean_length                 Mean fragment length of the samples 
                                       (if not accessible, use the reads' mean length)
         --std_length                  Std of the fragment length of the samples 
                                       (if not accessible, use the reads' length std or 10 if it is 0)

    Binning option:
        --bed                          Provides a .bed file with sites instead of processing it 
                                       from the reference genome with GATC_finder
        --fasta_bed                    Provide a fasta file of the bins 
        
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



params.skipFastp = false


include { fastp           } from "./nf_modules/fastp/main.nf"       
include { gatc_finder     } from "./nf_modules/gatc_finder/main.nf" addParams(gatc_finder_out: "${params.outdir}/sites_test/", gatc_finder: "--overlap")
include { multiqc         } from "./nf_modules/multiqc/main.nf"     addParams(multiqc_out: "${params.outdir}/reports/")
include { fasta_from_bed  } from "./nf_modules/bedtools/main.nf"    addParams(fasta_from_bed_out: "${params.outdir}/fasta/", fasta_from_bed: "-name")


params.single_end = false

params.mean_length = "150"
params.std_length = "10"


/* different inclusion of kallisto to pass different parameters if the reads are single end*/
if (params.single_end != false){

    if (params.mean_length == ""){
        error("In case of single end datas, please use --mean_length and --std_length to provide the mean and std of the fragments (if not known use the reads mean length and std)")
    } 
    if (params.std_length == "") {
        error("In case of single end datas, please use --mean_length and --std_length to provide the mean and std of the fragments (if not known use the reads mean length and std)")
    }

    include { index_fasta; mapping_fastq } from "./nf_modules/kallisto/main.nf" addParams(mapping_fastq_out: "${params.outdir}/counts/", mapping_fastq: "-l ${params.mean_length} -s ${params.std_length} --bias --bootstrap-samples 100")
} else {
    include { index_fasta; mapping_fastq } from "./nf_modules/kallisto/main.nf" addParams(mapping_fastq_out: "${params.outdir}/counts/")
}


params.fasta_bed = ""

params.bed = ""

/* creates channels according to the given arguments*/
if (params.fasta_bed == ""){

    

    if (params.bed != ""){
        channel
            .fromPath(params.bed)
            .ifEmpty { error "Cannot find any bed file matching: ${params.bed}" }
            .map { it -> [it.simpleName, it]}
            .set{bed}
    }
} else {
    channel
            .fromPath(params.fasta_bed)
            .ifEmpty { error "Cannot find any bed file matching: ${params.bed}" }
            .map { it -> [it.simpleName, it]}
            .set{fasta_bed}
}

params.fasta = ""

if (params.fasta_bed == ""){
    channel
        .fromPath(params.fasta)
        .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
        .map { it -> [it.simpleName, it]}
        .set {fasta_files}
}

params.fastq = ""

channel
    .fromFilePairs(params.fastq, size: -1)
    .set {fastq_files}



workflow {

    /*
    Classical wokrflow :
    - GATC_finder create a bed file containing all the GATC fragments of the given genome
    - Fasta_from_bed to convert the bed file into a fasta file
    - index_fasta to index the fasta file
    - Fastp to perform a qc on the reads
    - mapping_fastq to perform the pseumapping with kallisto
    - multiqc to gather the reports
    */


    if (params.fasta_bed == "") {

        /* skips GATC finder if --bed is given */ 
        if (params.bed != ""){

            fasta_from_bed(fasta_files,
                        bed)

        } else {

            gatc_finder(fasta_files)

            fasta_from_bed(fasta_files, 
                        gatc_finder.out.bed)
        }


        index_fasta(fasta_from_bed.out.fasta)

        /* skips fastp if skipsFastp is invoked*/
        if (params.skipFastp != false){


            mapping_fastq(index_fasta.out.index.collect(),
                        fastq_files)

            multiqc(mapping_fastq.out.report)

        } else {

            fastp(fastq_files)
            
            mapping_fastq(index_fasta.out.index.collect(),
                        fastp.out.fastq)

            report_mapping = mapping_fastq.out.report

            report_fastp = fastp.out.report

            multiqc(report_mapping.mix(report_fastp))
        }

    /* skips both GATC_finder and fasta_from_bed if --fasta_bed is given */
    } else {


        index_fasta(fasta_bed)


        if (params.skipFastp != false){


            mapping_fastq(index_fasta.out.index.collect(),
                        fastq_files)

            multiqc(mapping_fastq.out.report)

        } else {

            fastp(fastq_files)
            
            mapping_fastq(index_fasta.out.index.collect(),
                        fastp.out.fastq)

            report_mapping = mapping_fastq.out.report

            report_fastp = fastp.out.report

            multiqc(report_mapping.mix(report_fastp))


        }
    }
}