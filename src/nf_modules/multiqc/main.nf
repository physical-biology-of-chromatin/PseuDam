version = "1.9"
container_url = "lbmc/multiqc:${version}"

params.multiqc = ""
params.multiqc_out = ""
process multiqc {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  if (params.multiqc_out != "") {
    publishDir "results/${params.multiqc_out}", mode: 'copy'
  }

  input:
    path report 

  output:
    path "*multiqc_*", emit: report

  script:
"""
multiqc ${params.multiqc} -f .
"""
}
