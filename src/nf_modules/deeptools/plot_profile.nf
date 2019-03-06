params.matrix = "$baseDir/data/region_matrix/*.mat.gz"
params.title = "plot title"

log.info "matrix files : ${params.matrix}"
log.info "plot title : ${params.title}"

Channel
  .fromPath( params.matrix )
  .ifEmpty { error "Cannot find any matrix files matching: ${params.matrix}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { matrix_files }

process plot_profile {
  tag "$file_id"
  publishDir "results/mapping/region_matrix/", mode: 'copy'

  input:
    set file_id, file(matrix) from matrix_files

  output:
    set file_id, "*.pdf" into region_matrix

  script:
/*
see more option at
https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html
*/
"""
plotProfile -m ${matrix} \
  --plotFileFormat=pdf \
  -out ${file_id}.pdf \
  --plotType=fill \
  --perGroup \
  --plotTitle "${params.title}"
"""
}
