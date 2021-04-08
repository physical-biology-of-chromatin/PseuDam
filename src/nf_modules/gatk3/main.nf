version = "3.8.0"
container_url = "lbmc/gatk:${version}"

params.variant_calling = ""
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
"""
gatk3 -T HaplotypeCaller \
  -nct ${task.cpus} \
  ${params.variant_calling} \
  -R ${fasta} \
  -I ${bam} \
  -o ${file_id}.vcf
"""
}

params.filter_snp = ""
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
"""
gatk3 -T SelectVariants \
  -nct ${task.cpus} \
  ${params.filter_snp} \
  -R ${fasta} \
  -V ${vcf} \
  -selectType SNP \
  -o ${file_id}_snp.vcf
"""
}

params.filter_indels = ""
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
"""
gatk3 -T SelectVariants \
  -nct ${task.cpus} \
  ${params.filter_indels} \
  -R ${fasta} \
  -V ${vcf} \
  -selectType INDEL \
  -o ${file_id}_indel.vcf
"""
}

high_confidence_snp_filter = "(QD < 2.0) || (FS > 60.0) || (MQ < 40.0) || (MQRankSum < -12.5) || (ReadPosRankSum < -8.0) || (SOR > 4.0)"
params.high_confidence_snp = "--filterExpression \"${high_confidence_snp_filter}\" --filterName \"basic_snp_filter\""
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
"""
gatk3 -T VariantFiltration \
  -nct ${task.cpus} \
  -R ${fasta} \
  -V ${vcf} \
  ${params.high_confidence_snp} \
  -o ${file_id}_filtered_snp.vcf
"""
}

high_confidence_indel_filter = "QD < 3.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0"
params.high_confidence_indels = "--filterExpression \"${high_confidence_indel_filter}\" --filterName \"basic_indel_filter\""
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
"""
gatk3 -T VariantFiltration \
  -nct ${task.cpus} \
  -R ${fasta} \
  -V ${vcf} \
  ${params.high_confidence_indels} \
  -o ${file_id}_filtered_indel.vcf
"""
}

params.recalibrate_snp_table = ""
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
"""
gatk3 -T BaseRecalibrator \
  -nct ${task.cpus} \
  ${recalibrate_snp_table} \
  -R ${fasta} \
  -I ${bam} \
  -knownSites ${snp_file} \
  -knownSites ${indel_file} \
  -o recal_data_table
"""
}

params.recalibrate_snp = ""
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
"""
gatk3 -T PrintReads \
  --use_jdk_deflater \
  --use_jdk_inflater \
  ${recalibrate_snp} \
  -nct ${task.cpus} \
  -R ${fasta} \
  -I ${bam} \
  -BQSR recal_data_table \
  -o ${file_id}_recal.bam
"""
}

params.haplotype_caller = ""
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
"""
gatk3 -T HaplotypeCaller \
  -nct ${task.cpus} \
  ${params.haplotype_caller} \
  -R ${fasta} \
  -I ${bam} \
  -ERC GVCF \
  -variant_index_type LINEAR -variant_index_parameter 128000 \
  -o ${file_id}.gvcf
"""
}

params.gvcf_genotyping = ""
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
"""
gatk3 -T GenotypeGVCFs \
  -nct ${task.cpus} \
  ${params.gvcf_genotyping} \
  -R ${fasta} \
  -V ${gvcf} \
  -o ${file_id}_joint.vcf
"""
}

params.select_variants_snp = ""
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
"""
gatk3 -T SelectVariants \
  -nct ${task.cpus} \
  ${params.select_variants_snp} \
  -R ${fasta} \
  -V ${vcf} \
  -selectType SNP \
  -o ${file_id}_joint_snp.vcf
"""
}

params.select_variants_indels = ""
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
"""
gatk3 -T SelectVariants \
  -nct ${task.cpus} \
  ${params.select_variants_indels} \
  -R ${fasta} \
  -V ${vcf} \
  -selectType INDEL \
  -o ${file_id}_joint_indel.vcf
"""
}

params.personalized_genome = ""
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
  library = pick_library(file_id, library_list)
"""
gatk3 -T FastaAlternateReferenceMaker\
  ${params.personalized_genome} \
  -R ${reference} \
  -V ${vcf} \
  -o ${library}_genome.fasta
"""
}

