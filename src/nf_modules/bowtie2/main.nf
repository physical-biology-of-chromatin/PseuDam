version = "2.3.4.1"
container_url = "lbmc/bowtie2:${version}"

process index_fasta {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$fasta.baseName"

  input:
    path fasta

  output:
    path "*.index*", emit: index
    path "*_report.txt", emit: report

  script:
"""
bowtie2-build --threads ${task.cpus} \
  ${fasta} \
  ${fasta.baseName}.index &> \
  ${fasta.baseName}_bowtie2_index_report.txt

if grep -q "Error" ${fasta.baseName}_bowtie2_index_report.txt; then
  exit 1
fi
"""
}


process mapping_fastq {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$pair_id"

  input:
  path index
  tuple val(pair_id), path(reads)

  output:
  tuple val(pair_id), path("*.bam"), emit: bam
  path "*_report.txt", emit: report

  script:
  index_id = index[0]
  for (index_file in index) {
    if (index_file =~ /.*\.1\.bt2/ && !(index_file =~ /.*\.rev\.1\.bt2/)) {
        index_id = ( index_file =~ /(.*)\.1\.bt2/)[0][1]
    }
  }
if (reads instanceof List)
"""
bowtie2 --very-sensitive \
  -p ${task.cpus} \
  -x ${index_id} \
  -1 ${reads[0]} \
  -2 ${reads[1]} 2> \
  ${pair_id}_bowtie2_mapping_report_tmp.txt | \
  samtools view -Sb - > ${pair_id}.bam

if grep -q "Error" ${pair_id}_bowtie2_mapping_report_tmp.txt; then
  exit 1
fi
tail -n 19 ${pair_id}_bowtie2_mapping_report_tmp.txt > \
  ${pair_id}_bowtie2_mapping_report.txt
"""
else
"""
bowtie2 --very-sensitive \
  -p ${task.cpus} \
  -x ${index_id} \
  -U ${reads} 2> \
  ${reads.baseName}_bowtie2_mapping_report_tmp.txt | \
  samtools view -Sb - > ${reads.baseName}.bam

if grep -q "Error" ${reads.baseName}_bowtie2_mapping_report_tmp.txt; then
  exit 1
fi
tail -n 19 ${reads.baseName}_bowtie2_mapping_report_tmp.txt > \
  ${reads.baseName}_bowtie2_mapping_report.txt
"""
}

process mapping_fastq_pairedend {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$pair_id"

  input:
  path index
  tuple val(pair_id), path(reads)

  output:
  tuple val(pair_id), path("*.bam"), emit: bam
  path "*_report.txt", emit: report

  script:
  index_id = index[0]
  for (index_file in index) {
    if (index_file =~ /.*\.1\.bt2/ && !(index_file =~ /.*\.rev\.1\.bt2/)) {
        index_id = ( index_file =~ /(.*)\.1\.bt2/)[0][1]
    }
  }
"""
bowtie2 --very-sensitive \
  -p ${task.cpus} \
  -x ${index_id} \
  -1 ${reads[0]} \
  -2 ${reads[1]} 2> \
  ${pair_id}_bowtie2_mapping_report_tmp.txt | \
  samtools view -Sb - > ${pair_id}.bam

if grep -q "Error" ${pair_id}_bowtie2_mapping_report_tmp.txt; then
  exit 1
fi
tail -n 19 ${pair_id}_bowtie2_mapping_report_tmp.txt > \
  ${pair_id}_bowtie2_mapping_report.txt
"""
}


process mapping_fastq_singleend {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
  path index
  tuple val(file_id), path(reads)

  output:
  tuple val(file_id), path("*.bam"), emit: bam
  path "*_report.txt", emit: report

  script:
  index_id = index[0]
  for (index_file in index) {
    if (index_file =~ /.*\.1\.bt2/ && !(index_file =~ /.*\.rev\.1\.bt2/)) {
        index_id = ( index_file =~ /(.*)\.1\.bt2/)[0][1]
    }
  }
"""
bowtie2 --very-sensitive \
  -p ${task.cpus} \
  -x ${index_id} \
  -U ${reads} 2> \
  ${reads.baseName}_bowtie2_mapping_report_tmp.txt | \
  samtools view -Sb - > ${reads.baseName}.bam

if grep -q "Error" ${reads.baseName}_bowtie2_mapping_report_tmp.txt; then
  exit 1
fi
tail -n 19 ${reads.baseName}_bowtie2_mapping_report_tmp.txt > \
  ${reads.baseName}_bowtie2_mapping_report.txt
"""
}
