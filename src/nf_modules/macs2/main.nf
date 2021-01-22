version = "2.1.2"
container_url = "lbmc/macs2:${version}"

process peak_calling {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "${file_id}"
  publishDir "results/peak_calling/${file_id}", mode: 'copy'

  input:
    tuple val(file_id), path(bam_ip), path(bam_control)

  output:
    path "*", emit: peak
    path "*_report.txt", emit: report

  script:
/* remove --nomodel option for real dataset */
"""
macs2 callpeak \
  --treatment ${file_ip} \
  --call-summits "True"\
  --control ${file_control} \
  --keep-dup "auto" \
  --name ${file_id} \
  --gsize ${params.genome_size} 2> \
  ${file_ip}_macs2_report.txt

if grep -q "ERROR" ${file_ip}_macs2_report.txt; then
  echo "MACS2 error"
  exit 1
fi
"""
}

