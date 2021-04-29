version = "dd69682"
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
    tuple val(file_id), path(bam)
    tuple val(gtf_id), path(gtf)

  output:
    tuple val(file_id), path("${bam.simpleName}.bin"), emit: bin
    tuple val(gtf_id), path("${gtf.simpleName}_transcripts_lengths.tsv"), emit: tsv

  script:
"""
awk -F"\t" '
\$3=="exon" 
    {
        ID=substr(\$9, length(\$9)-16, 15); 
        L[ID]+=\$5-\$4+1
    } 
END{
    for(i in L)
        {print i"\t"L[i]}
    }
' ${gtf} ${gtf.simpleName}_transcripts_lengths.tsv
alntools bam2ec ${params.bam2sec} -t ${gtf.simpleName}_transcripts_lengths.tsv -c 8 ${bam} ${bam.simpleName}.bin
"""
}