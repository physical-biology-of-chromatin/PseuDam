version = "2.25.0"
container_url = "lbmc/bedtools:${version}"

params.fasta_from_bed = "-name"
process fasta_from_bed {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${bed.baseName}"

  input:
  path fasta
  path bed

  output:
  tuple val(bed.baseName), path("*_extracted.fasta"), emit: fasta

  script:
"""
bedtools getfasta ${params.fasta_from_bed} \
-fi ${fasta} -bed ${bed} -fo ${bed.baseName}_extracted.fasta
"""
}

params.merge_bed = ""
process merge_bed {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${bed.baseName}"

  input:
  path bed

  output:
  tuple val(bed[0].simpleName), path("*_merged.fasta"), emit: bed

  script:
"""
bedtools merge ${params.merge_bed} -i ${bed} > ${bed[0].simpleName}_merged.bed
"""
}

params.bam_to_fastq_singleend = ""
process bam_to_fastq_singleend {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${bam_id}"

  input:
  tuple val(bam_id), path(bam)

  output:
  tuple val(bam_id), path("*.fastq"), emit: fastq

  script:
"""
bedtools bamtofastq \
  ${params.bam_to_fastq_singleend} \
  -i ${bam} -fq ${bam.baseName}.fastq
"""
}

params.bam_to_fastq_pairedend = ""
process bam_to_fastq_pairedend {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${bam_id}"

  input:
  tuple val(bam_id), path(bam)

  output:
  tuple val(bam_id), path("*.fastq"), emit: fastq

  script:
"""
bedtools bamtofastq \
  ${params.bam_to_fastq_pairedend} \
  -i ${bam} -fq ${bam.baseName}_R1.fastq -fq2 ${bam.baseName}_R2.fastq
"""
}

params.bam_to_bedgraph = ""
process bam_to_bedgraph {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${bam_id}"

  input:
  tuple val(bam_id), path(bam)

  output:
  tuple val(bam_id), path("*.bg"), emit: bedgraph

  script:
"""
bedtools genomecov \
  ${params.bam_to_bedgraph} \
  -ibam ${bam} \
  -bg > ${bam.simpleName}.bg
"""
}
