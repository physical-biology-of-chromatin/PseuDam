params.bw = "$baseDir/data/bigwig/*.bw"
params.bed = "$baseDir/data/annot/*.bed"

log.info "bigwig files : ${params.bw}"
log.info "bed files : ${params.bed}"

Channel
  .fromPath( params.bw )
  .ifEmpty { error "Cannot find any bigwig files matching: ${params.bw}" }
  .set { bw_files }

Channel
  .fromPath( params.bed )
  .ifEmpty { error "Cannot find any bed files matching: ${params.bed}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { bed_files }

process compute_matrix {
  tag "$bed_file_id"
  cpus 4
  publishDir "results/mapping/region_matrix/", mode: 'copy'

  input:
    file bw from bw_files.collect()
    set bed_file_id, file(bed) from bed_files.collect()

  output:
    set bed_file_id, "*.mat.gz" into region_matrix

  script:
"""
computeMatrix scale-regions -S ${bw} \
  -p ${task.cpus} \
  -R ${bed} \
  --beforeRegionStartLength 100 \
  --afterRegionStartLength 100 \
  -o ${bed_file_id}.mat.gz
"""
}
