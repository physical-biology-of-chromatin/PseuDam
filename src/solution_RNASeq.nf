nextflow.enable.dsl=2

/*
./nextflow src/solution_RNASeq.nf --fastq "data/tiny_dataset/fastq/tiny2_R{1,2}.fastq.gz" --fasta "data/tiny_dataset/fasta/tiny_v2_10.fasta" --bed "data/tiny_dataset/annot/tiny.bed" -profile docker
*/

log.info "fastq files : ${params.fastq}"
log.info "fasta file : ${params.fasta}"
log.info "bed file : ${params.bed}"

Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
  .set { fasta_files }
Channel
  .fromPath( params.bed )
  .ifEmpty { error "Cannot find any bed files matching: ${params.bed}" }
  .set { bed_files }
Channel
  .fromFilePairs( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

include { adaptor_removal } from './nf_modules/cutadapt/main'
include { trimming } from './nf_modules/urqt/main'
include { fasta_from_bed } from './nf_modules/bedtools/main'
include { index_fasta; mapping_fastq } from './nf_modules/kallisto/main'

workflow {
    adaptor_removal(fastq_files)
    trimming(adaptor_removal.out.fastq)
    fasta_from_bed(fasta_files, bed_files)
    index_fasta(fasta_from_bed.out.fasta)
    mapping_fastq(index_fasta.out.index.collect(), trimming.out.fastq)
}

