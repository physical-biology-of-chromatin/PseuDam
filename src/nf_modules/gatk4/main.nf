version = "4.2.0.0"
container_url = "broadinstitute/gatk:${version}"

process variant_calling {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam), path(bai)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.vcf"), emit: vcf

  script:
  xmx_memory = task.memory - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" -T HaplotypeCaller \
  -nct ${task.cpus} \
  -R ${fasta} \
  -I ${bam} \
  -o ${file_id}.vcf
"""
}

process filter_snp {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_snp.vcf"), emit: vcf
  script:
  xmx_memory = task.memory - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" -T SelectVariants \
  -nct ${task.cpus} \
  -R ${fasta} \
  -V ${vcf} \
  -selectType SNP \
  -o ${file_id}_snp.vcf
"""
}

process filter_indels {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_indel.vcf"), emit: vcf
  script:
  xmx_memory = task.memory - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G"-T SelectVariants \
  -nct ${task.cpus} \
  -R ${fasta} \
  -V ${vcf} \
  -selectType INDEL \
  -o ${file_id}_indel.vcf
"""
}

high_confidence_snp_filter = "(QD < 2.0) || (FS > 60.0) || (MQ < 40.0) || (MQRankSum < -12.5) || (ReadPosRankSum < -8.0) || (SOR > 4.0)"

process high_confidence_snp {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_snp.vcf"), emit: vcf
  script:
  xmx_memory = task.memory - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G"-T VariantFiltration \
  -nct ${task.cpus} \
  -R ${fasta} \
  -V ${vcf} \
  --filterExpression "${high_confidence_snp_filter}" \
  --filterName "basic_snp_filter" \
  -o ${file_id}_filtered_snp.vcf
"""
}

high_confidence_indel_filter = "QD < 3.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0"

process high_confidence_indels {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_indel.vcf"), emit: vcf
  script:
  xmx_memory = task.memory - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" -T VariantFiltration \
  -nct ${task.cpus} \
  -R ${fasta} \
  -V ${vcf} \
  --filterExpression "${high_confidence_indel_filter}" \
  --filterName "basic_indel_filter" \
  -o ${file_id}_filtered_indel.vcf
"""
}

process recalibrate_snp_table {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(snp_file), path(indel_file), path(bam), path(bam_idx)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("recal_data_table"), emit: recal_table
  script:
  xmx_memory = task.memory - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" -T BaseRecalibrator \
  -nct ${task.cpus} \
  -R ${fasta} \
  -I ${bam} \
  -knownSites ${snp_file} \
  -knownSites ${indel_file} \
  -o recal_data_table
"""
}

process recalibrate_snp {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(snp_file), path(indel_file), path(bam), path(bam_idx)
    tuple val(table_id), path(recal_data_table)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.bam"), emit: bam
  script:
  xmx_memory = task.memory - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" -T PrintReads \
  --use_jdk_deflater \
  --use_jdk_inflater \
  -nct ${task.cpus} \
  -R ${fasta} \
  -I ${bam} \
  -BQSR recal_data_table \
  -o ${file_id}_recal.bam
"""
}

process haplotype_caller {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.gvcf"), emit: gvcf
  script:
  xmx_memory = task.memory - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" -T HaplotypeCaller \
  -nct ${task.cpus} \
  -R ${fasta} \
  -I ${bam} \
  -ERC GVCF \
  -variant_index_type LINEAR -variant_index_parameter 128000 \
  -o ${file_id}.gvcf
"""
}

process gvcf_genotyping {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(gvcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.vcf"), emit: vcf
  script:
  xmx_memory = task.memory - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" -T GenotypeGVCFs \
  -nct ${task.cpus} \
  -R ${fasta} \
  -V ${gvcf} \
  -o ${file_id}_joint.vcf
"""
}

process select_variants_snp {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_joint_snp.vcf"), emit: vcf
  script:
  xmx_memory = task.memory - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}GG" -T SelectVariants \
  -nct ${task.cpus} \
  -R ${fasta} \
  -V ${vcf} \
  -selectType SNP \
  -o ${file_id}_joint_snp.vcf
"""
}

process select_variants_indels {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_joint_indel.vcf"), emit: vcf
  script:
  xmx_memory = task.memory - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" -T SelectVariants \
  -nct ${task.cpus} \
  -R ${fasta} \
  -V ${vcf} \
  -selectType INDEL \
  -o ${file_id}_joint_indel.vcf
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
  xmx_memory = task.memory - ~/\s*GB/
"""
gatk --java-options "-Xmx${xmx_memory}G" -T FastaAlternateReferenceMaker\
  -R ${reference} \
  -V ${vcf} \
  -o ${file_id[0]}_genome.fasta
"""
}

