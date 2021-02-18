version = "3.0.0a6"
container_url = "lbmc/macs3:${version}"

macs3_genome_size=3e9

process peak_calling {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "${file_id}"

  input:
    tuple val(file_id), path(bam_ip), path(bam_control)

  output:
    path "*", emit: peak
    path "*_report.txt", emit: report

  script:
/* remove --nomodel option for real dataset */
"""
macs2 callpeak \
  --treatment ${bam_ip} \
  --call-summits "True"\
  --control ${bam_control} \
  --keep-dup "auto" \
  --name ${bam_ip.simpleName} \
  --gsize ${macs3_genome_size} 2> \
  ${bam_ip.simpleName}_macs3_report.txt

if grep -q "ERROR" ${bam_ip.simpleName}_macs3_report.txt; then
  echo "MACS3 error"
  exit 1
fi
"""
}

process peak_calling_bg {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "${file_id}"

  input:
    tuple val(file_id), path(bg_ip), path(bg_control)

  output:
    path "*", emit: peak
    path "*_report.txt", emit: report

  script:
/* remove --nomodel option for real dataset */
"""
awk '{print \$1"\t"\$2"\t"\$3"\t.\t+\t"\$4}' ${bg_ip} > \
  ${bg_ip.simpleName}.bed
awk '{print \$1"\t"\$2"\t"\$3"\t.\t+\t"\$4}' ${bg_control} > \
  ${bg_control.simpleName}.bed
macs2 callpeak \
  --treatment ${bg_ip.simpleName}.bed \
  --call-summits "True"\
  --control ${bg_control.simpleName}.bed \
  --keep-dup "auto" \
  --name ${bg_ip.simpleName} \
  --gsize ${macs3_genome_size} 2> \
  ${bg_ip.simpleName}_macs3_report.txt

if grep -q "ERROR" ${bg_ip.simpleName}_macs3_report.txt; then
  echo "MACS3 error"
  exit 1
fi
"""
}

