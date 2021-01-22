version = "1.2.2"
container_url = "lbmc/bowtie:${version}"

process index_fasta {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$fasta.baseName"
  publishDir "results/mapping/index/", mode: 'copy'

  input:
    path fasta

  output:
    path "*.index*", emit: index
    path "*_report.txt", emit: report

  script:
"""
bowtie-build --threads ${task.cpus} \
  -f ${fasta} ${fasta.baseName}.index &> \
  ${fasta.baseName}_bowtie_index_report.txt

if grep -q "Error" ${fasta.baseName}_bowtie_index_report.txt; then
  exit 1
fi
"""
}


process mapping_fastq_pairedend {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$pair_id"
  publishDir "results/mapping/bams/", mode: 'copy'

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
# -v specify the max number of missmatch, -k the number of match reported per
# reads
bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_id} \
  -1 ${reads[0]} -2 ${reads[1]} 2> \
  ${pair_id}_bowtie_report_tmp.txt | \
  samtools view -Sb - > ${pair_id}.bam

if grep -q "Error" ${pair_id}_bowtie_report_tmp.txt; then
  exit 1
fi
tail -n 19 ${pair_id}_bowtie_report_tmp.txt > \
  ${pair_id}_bowtie_mapping_report.txt
"""
}


process mapping_fastq_singleend {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  publishDir "results/mapping/bams/", mode: 'copy'

  input:
  path index
  tuple val(file_id), path(reads)

  output:
  set file_id, "*.bam", emit: bam
  file "*_report.txt", emit: report

  script:
  index_id = index[0]
  for (index_file in index) {
    if (index_file =~ /.*\.1\.bt2/ && !(index_file =~ /.*\.rev\.1\.bt2/)) {
        index_id = ( index_file =~ /(.*)\.1\.bt2/)[0][1]
    }
  }
"""
bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_id} \
  -q ${reads} 2> \
  ${file_id}_bowtie_report_tmp.txt | \
  samtools view -Sb - > ${file_id}.bam

if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
  exit 1
fi
tail -n 19 ${file_id}_bowtie_report_tmp.txt > \
  ${file_id}_bowtie_mapping_report.txt
"""
}
