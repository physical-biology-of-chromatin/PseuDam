version = "4.2.0.0"
container_url = "broadinstitute/gatk:${version}"

process variant_calling {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam), path(bai)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.vcf"), emit: vcf

  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" HaplotypeCaller \
  -R ${fasta} \
  -I ${bam} \
  -O ${file_id}.vcf
"""
}

process filter_snp {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_snp.vcf"), emit: vcf
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" SelectVariants \
  -R ${fasta} \
  -V ${vcf} \
  -select-type SNP \
  -O ${file_id}_snp.vcf
"""
}

process filter_indels {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_indel.vcf"), emit: vcf
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" SelectVariants \
  -R ${fasta} \
  -V ${vcf} \
  -select-type INDEL \
  -O ${file_id}_indel.vcf
"""
}

high_confidence_snp_filter = "(QD < 2.0) || (FS > 60.0) || (MQ < 40.0) || (MQRankSum < -12.5) || (ReadPosRankSum < -8.0) || (SOR > 4.0)"

process high_confidence_snp {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_snp.vcf"), emit: vcf
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" VariantFiltration \
  -R ${fasta} \
  -V ${vcf} \
  --filter-expression "${high_confidence_snp_filter}" \
  --filter-name "basic_snp_filter" \
  -O ${file_id}_filtered_snp.vcf
"""
}

high_confidence_indel_filter = "QD < 3.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0"

process high_confidence_indels {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_indel.vcf"), emit: vcf
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" VariantFiltration \
  -R ${fasta} \
  -V ${vcf} \
  --filter-expression "${high_confidence_indel_filter}" \
  --filter-name "basic_indel_filter" \
  -O ${file_id}_filtered_indel.vcf
"""
}

process recalibrate_snp_table {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(snp_file), path(indel_file), path(bam), path(bam_idx)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("recal_data_table"), emit: recal_table
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" IndexFeatureFile \
  -I ${snp_file}
gatk --java-options "-Xmx${xmx_memory}G" IndexFeatureFile \
  -I ${indel_file}
gatk --java-options "-Xmx${xmx_memory}G" BaseRecalibrator \
  -R ${fasta} \
  -I ${bam} \
  -known-sites ${snp_file} \
  -known-sites ${indel_file} \
  -O recal_data_table
"""
}

process recalibrate_snp {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(snp_file), path(indel_file), path(bam), path(bam_idx), path(recal_table)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.bam"), emit: bam
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" ApplyBQSR \
  -R ${fasta} \
  -I ${bam} \
  --bqsr-recal-file recal_data_table \
  -O ${file_id}_recal.bam
"""
}

process haplotype_caller {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.gvcf"), emit: gvcf
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" HaplotypeCaller \
  -R ${fasta} \
  -I ${bam} \
  -ERC GVCF \
  -O ${file_id}.gvcf
"""
}

process gvcf_genotyping {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(gvcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.vcf"), emit: vcf
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" GenotypeGVCFs \
  -R ${fasta} \
  -V ${gvcf} \
  -O ${file_id}_joint.vcf
"""
}

process select_variants_snp {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_joint_snp.vcf"), emit: vcf
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}GG" SelectVariants \
  -R ${fasta} \
  -V ${vcf} \
  -select-type SNP \
  -O ${file_id}_joint_snp.vcf
"""
}

process select_variants_indels {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_joint_indel.vcf"), emit: vcf
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" SelectVariants \
  -R ${fasta} \
  -V ${vcf} \
  -select-type INDEL \
  -O ${file_id}_joint_indel.vcf
"""
}

process personalized_genome {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_genome.fasta"), emit: fasta

  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" FastaAlternateReferenceMaker\
  -R ${reference} \
  -V ${vcf} \
  -O ${file_id[0]}_genome.fasta
"""
}

