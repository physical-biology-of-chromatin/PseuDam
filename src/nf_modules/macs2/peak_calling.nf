params.genome_size = "hs"
params.control_tag = "control"
log.info "bam files : ${params.bam}"
log.info "genome size : ${params.genome_size}"

Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any bam files matching: ${params.bam}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { bam_files }

/* split bam Channel into control and ip if "control" keyword is detected*/
bam_files_control = Channel.create()
bam_files_ip = Channel.create()
bam_files.choice(
  bam_files_control,
  bam_files_ip ) { a -> a[0] =~ /.*${params.control_tag}.*/ ? 0 : 1 }

process peak_calling {
  tag "${file_id}"
  publishDir "results/peak_calling/${file_id}", mode: 'copy'

  input:
    set file_id, file(file_ip) from bam_files_ip
    set file_id_control, file(file_control) from bam_files_control.collect()

  output:
    file "*" into peak_output
    file "*_report.txt" into peak_calling_report

  script:
/* remove --nomodel option for real dataset */
"""
macs2 callpeak \
  --nomodel \
  --treatment ${file_ip} \
  --control ${file_control} \
  --name ${file_id} \
  --gsize ${params.genome_size} 2> \
${file_ip}_macs2_report.txt

if grep -q "ERROR" ${file_ip}_macs2_report.txt; then
  echo "MACS2 error"
  exit 1
fi
"""
}
