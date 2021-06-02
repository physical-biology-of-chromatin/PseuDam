nextflow.enable.dsl=2

/*
Testing pipeline for marseq scRNASeq analysis
*/

include {
  mapping;
} from "./nf_modules/bwa/main.nf"

include {
  variant_calling_out;
} from "./nf_modules/gatk4/main.nf" addParams(
  variant_calling_out: "vcf/",
)

channel
  .fromFilePairs( params.fastq, size: -1)
  .set { fastq_files }
channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
  .map { it -> [it.simpleName, it]}
  .set { fasta_files }

workflow {
  mapping(fasta, fastq_files)
  germline_cohort_data_variant_calling(mapping.out.bam, fasta)
}
