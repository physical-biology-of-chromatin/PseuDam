version = "0.8.0"
container_url = "lbmc/agat:${version}"

params.gff_to_bed = ""
params.gff_to_bed_out = ""
process gff_to_bed {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.gff_to_bed_out != "") {
    publishDir "results/${params.gff_to_bed_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(gff)
  output:
    tuple val(file_id), path("*.bed"), emit: bed

  script:
"""
agat_convert_sp_gff2bed.pl --gff ${gff} -o ${gff.simpleName}.bed
"""
}