version = "dd96682"
container_url = "lbmc/alntools:${version}"

params.bam2ec = ""
params.bam2ec_out = ""
process bam2ec {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.bam2ec_out != "") {
    publishDir "results/${params.bam2ec_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam), path(bam_idx)
    tuple val(transcripts_lengths_id), path(transcripts_lengths)

  output:
    tuple val(file_id), path("${bam.simpleName}.bin"), emit: bin
    tuple val(gtf_id), path("${transcripts_lengths}"), emit: tsv

  script:
"""
cp ${bam} file_bam.bam
cp ${bam_idx} file_bam.bam.bai
alntools bam2ec ${params.bam2ec} -t ${transcripts_lengths} file_bam.bam ${bam.simpleName}.bin
"""
}

params.gtf_to_transcripts_lengths = ""
params.gtf_to_transcripts_lengths_out = ""
process gtf_to_transcripts_lengths {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.gtf_to_transcripts_lengths != "") {
    publishDir "results/${params.gtf_to_transcripts_lengths}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(gtf)

  output:
    tuple val(file_id), path("${gtf.simpleName}_transcripts_lengths.tsv"), emit: tsv

  script:
"""
awk -F"[\\t;]" '
\$3=="exon" {
        ID=gensub(/transcript_id \\"(.*)\\"/, "\\\\1", "g", \$11); 
        LEN[ID]+=\$5-\$4+1;
    } 
END{
    for(i in LEN)
        {print i"\\t"LEN[i]}
    }
' ${gtf} > ${gtf.simpleName}_transcripts_lengths.tsv
"""
}
