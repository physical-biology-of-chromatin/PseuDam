version = "2.25.0"
container_url = "lbmc/bedtools:${version}"

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
bedtools getfasta -name \
-fi ${fasta} -bed ${bed} -fo ${bed.baseName}_extracted.fasta
"""
}

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
bedtools merge -i ${bed} > ${bed[0].simpleName}_merged.bed
"""
}

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
-i ${bam} -fq ${bam.baseName}.fastq
"""
}

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
-i ${bam} -fq ${bam.baseName}_R1.fastq -fq2 ${bam.baseName}_R2.fastq
"""
}

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
  -ibam ${bam} \
  -bg > ${bam.simpleName}.bg
"""
}
