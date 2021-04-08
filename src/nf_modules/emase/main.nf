version = "0.10.16"
container_url = "lbmc/emase:${version}"

params.personalised_transcriptome = ""

process personalised_transcriptome {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(fasta)
    tuple val(gtf_id), path(gtf)

  output:
    tuple val(file_id), path("${fasta.simpleName}.*"), emit: index
    tuple val(file_id), path("*_bwa_report.txt"), emit: report

  script:
"""
prepare-emase ${personalised_transcriptome} -G ${REF_FASTA} -g ${REF_GTF} -o ${REF_DIR} -m --no-bowtie-index
// ${REF_DIR}/emase.transcriptome.fa
// ${REF_DIR}/emase.transcriptome.info
// ${REF_DIR}/emase.gene2transcripts.tsv
prepare-emase -G ${SAMPLE_DIR}/L.fa,${SAMPLE_DIR}/R.fa -s L,R -o ${SAMPLE_DIR}
"""
}