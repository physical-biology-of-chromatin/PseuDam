nextflow.enable.dsl=2

/*
Testing pipeline for marseq scRNASeq analysis
*/

include { index_fasta; count } from "./nf_modules/kb/main.nf" addParams(
  kb_protocol: "marsseq",
  count_out: "quantification/"
)

params.fasta = "http://ftp.ensembl.org/pub/release-94/fasta/gallus_gallus/dna/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa.gz"
params.cdna = "http://ftp.ensembl.org/pub/release-94/fasta/gallus_gallus/cdna/Gallus_gallus.Gallus_gallus-5.0.cdna.all.fa.gz"
params.fastq = "data/CF42_45/*/*R{1,2}.fastq.gz"
params.gtf = "http://ftp.ensembl.org/pub/release-94/gtf/gallus_gallus/Gallus_gallus.Gallus_gallus-5.0.94.gtf.gz"
params.transcript_to_gene = ""
params.whitelist = "data/expected_whitelist.txt"
params.config = "data/marseq_flexi_splitter.yaml"

log.info "fastq files: ${params.fastq}"
log.info "fasta file : ${params.fasta}"
log.info "gtf file : ${params.gtf}"
log.info "transcript_to_gene file : ${params.transcript_to_gene}"
log.info "whitelist file : ${params.whitelist}"

channel
  .fromFilePairs( params.fastq, size: -1)
  .set { fastq_files }
channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
  .map { it -> [it.simpleName, it]}
  .set { fasta_files }
channel
  .fromPath( params.cdna )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.cdna}" }
  .map { it -> [it.simpleName, it]}
  .set { cdna_files }
channel
  .fromPath( params.gtf )
  .ifEmpty { error "Cannot find any gtf files matching: ${params.gtf}" }
  .map { it -> [it.simpleName, it]}
  .set { gtf_files }
if (params.whitelist == "") {
  channel.empty()
    .set { whitelist_files }
} else {
  channel
    .fromPath( params.whitelist )
    .map { it -> [it.simpleName, it]}
    .set { whitelist_files }
}
channel
  .fromPath( params.config )
  .ifEmpty { error "Cannot find any config files matching: ${params.config}" }
  .map { it -> [it.simpleName, it]}
  .set { config_files }

workflow {
  index_fasta(fasta_files, cdna_files, gtf_files)
  count(index_fasta.out.index, fastq_files, index_fasta.out.t2g, whitelist_files, config_files)
}
