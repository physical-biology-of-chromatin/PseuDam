version = "1.9"
container_url = "lbmc/multiqc:${version}"

params.multiqc = ""
process multiqc {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  publishDir "results/QC/", mode: 'copy'

  input:
    path report

  output:
    path "*multiqc_*", emit: report

  script:
"""
multiqc ${params.multiqc} -f .
"""
}
