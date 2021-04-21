version = "4.2.0.0"
container_url = "broadinstitute/gatk:${version}"

params.variant_calling = ""
params.variant_calling_out = ""
process variant_calling {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.variant_calling_out != "") {
    publishDir "results/${params.variant_calling_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam), path(bai)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.vcf"), emit: vcf

  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk --java-options "-Xmx${xmx_memory}G" HaplotypeCaller \
  ${params.variant_calling} \
  -R ${fasta} \
  -I ${bam} \
  -O ${bam.simpleName}.vcf
"""
}

params.filter_snp = ""
params.filter_snp_out = ""
process filter_snp {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.filter_snp_out != "") {
    publishDir "results/${params.filter_snp_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_snp.vcf"), emit: vcf
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk --java-options "-Xmx${xmx_memory}G" SelectVariants \
  ${params.filter_snp} \
  -R ${fasta} \
  -V ${vcf} \
  -select-type SNP \
  -O ${vcf.simpleName}_snp.vcf
"""
}

params.filter_indels = ""
params.filter_indels_out = ""
process filter_indels {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.filter_indels_out != "") {
    publishDir "results/${params.filter_indels_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_indel.vcf"), emit: vcf
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk --java-options "-Xmx${xmx_memory}G" SelectVariants \
  ${params.filter_indels} \
  -R ${fasta} \
  -V ${vcf} \
  -select-type INDEL \
  -O ${vcf.simpleName}_indel.vcf
"""
}

params.high_confidence_snp_filter = "(QD < 2.0) || (FS > 60.0) || (MQ < 40.0) || (MQRankSum < -12.5) || (ReadPosRankSum < -8.0) || (SOR > 4.0)"
params.high_confidence_snp = "--filter-expression \"${params.high_confidence_snp_filter}\" --filter-name \"basic_snp_filter\""
params.high_confidence_snp_out = ""
process high_confidence_snp {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.high_confidence_snp_out != "") {
    publishDir "results/${params.high_confidence_snp_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_snp.vcf"), emit: vcf
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk --java-options "-Xmx${xmx_memory}G" VariantFiltration \
  -R ${fasta} \
  -V ${vcf} \
  ${params.high_confidence_snp} \
  -O ${vcf.simpleName}_filtered_snp.vcf
"""
}

params.high_confidence_indel_filter = "QD < 3.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0"
params.high_confidence_indels = "--filter-expression \"${params.high_confidence_indel_filter}\" --filter-name \"basic_indel_filter\""
params.high_confidence_indels_out = ""
process high_confidence_indels {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.high_confidence_indels_out != "") {
    publishDir "results/${params.high_confidence_indels_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_indel.vcf"), emit: vcf
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk --java-options "-Xmx${xmx_memory}G" VariantFiltration \
  -R ${fasta} \
  -V ${vcf} \
  ${params.high_confidence_indels} \
  -O ${vcf.simpleName}_filtered_indel.vcf
"""
}

params.recalibrate_snp_table = ""
params.recalibrate_snp_table_out = ""
process recalibrate_snp_table {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.recalibrate_snp_table_out != "") {
    publishDir "results/${params.recalibrate_snp_table_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(snp_file), path(indel_file), path(bam), path(bam_idx)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("recal_data_table"), emit: recal_table
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk --java-options "-Xmx${xmx_memory}G" IndexFeatureFile \
  -I ${snp_file}
gatk --java-options "-Xmx${xmx_memory}G" IndexFeatureFile \
  -I ${indel_file}
gatk --java-options "-Xmx${xmx_memory}G" BaseRecalibrator \
  ${params.recalibrate_snp_table} \
  -R ${fasta} \
  -I ${bam} \
  -known-sites ${snp_file} \
  -known-sites ${indel_file} \
  -O recal_data_table
"""
}

params.recalibrate_snp = ""
params.recalibrate_snp_out = ""
process recalibrate_snp {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.recalibrate_snp_out != "") {
    publishDir "results/${params.recalibrate_snp_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(snp_file), path(indel_file), path(bam), path(bam_idx), path(recal_table)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.bam"), emit: bam
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk --java-options "-Xmx${xmx_memory}G" ApplyBQSR \
  ${params.recalibrate_snp} \
  -R ${fasta} \
  -I ${bam} \
  --bqsr-recal-file recal_data_table \
  -O ${bam.simpleName}_recal.bam
"""
}

params.haplotype_caller = ""
params.haplotype_caller_out = ""
process haplotype_caller {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.haplotype_caller_out != "") {
    publishDir "results/${params.haplotype_caller_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.gvcf"), emit: gvcf
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk --java-options "-Xmx${xmx_memory}G" HaplotypeCaller \
  ${params.haplotype_caller} \
  -R ${fasta} \
  -I ${bam} \
  -ERC GVCF \
  -O ${bam.simpleName}.gvcf
"""
}

params.gvcf_genotyping = ""
params.gvcf_genotyping_out = ""
process gvcf_genotyping {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.gvcf_genotyping_out != "") {
    publishDir "results/${params.gvcf_genotyping_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(gvcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.vcf.gz"), emit: vcf
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk --java-options "-Xmx${xmx_memory}G" GenotypeGVCFs \
  ${params.gvcf_genotyping} \
  -R ${fasta} \
  -V ${gvcf} \
  -O ${gvcf.simpleName}_joint.vcf.gz
"""
}

params.select_variants_snp = ""
params.select_variants_snp_out = ""
process select_variants_snp {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.select_variants_snp_out != "") {
    publishDir "results/${params.select_variants_snp_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_joint_snp.vcf"), emit: vcf
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk --java-options "-Xmx${xmx_memory}GG" SelectVariants \
  ${params.select_variants_snp} \
  -R ${fasta} \
  -V ${vcf} \
  -select-type SNP \
  -O ${vcf.simpleName}_joint_snp.vcf
"""
}

params.select_variants_indels = ""
params.select_variants_indels_out = ""
process select_variants_indels {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.select_variants_indels_out != "") {
    publishDir "results/${params.select_variants_indels_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_joint_indel.vcf"), emit: vcf
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk --java-options "-Xmx${xmx_memory}G" SelectVariants \
  ${params.select_variants_indels} \
  -R ${fasta} \
  -V ${vcf} \
  -select-type INDEL \
  -O ${file_prefix}_joint_indel.vcf
"""
}

params.personalized_genome = ""
params.personalized_genome_out = ""
process personalized_genome {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.personalized_genome_out != "") {
    publishDir "results/${params.personalized_genome_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_genome.fasta"), emit: fasta

  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk --java-options "-Xmx${xmx_memory}G" FastaAlternateReferenceMaker\
  ${params.personalized_genome} \
  -R ${reference} \
  -V ${vcf} \
  -O ${vcf.simpleName}_genome.fasta
"""
}



