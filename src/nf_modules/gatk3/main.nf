version = "3.8.0"
container_url = "lbmc/gatk:${version}"

process variant_calling {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam), path(bai)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), "*.vcf", emit: vcf

  script:
"""
gatk3 -T HaplotypeCaller \
  -nct ${task.cpus} \
  -R ${fasta} \
  -I ${bam} \
  -o ${file_id}.vcf
"""
}

process filter_snp {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(variants)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_snp.vcf"), emit: vcf
  script:
"""
gatk3 -T SelectVariants \
  -nct ${task.cpus} \
  -R ${fasta} \
  -V ${variants} \
  -selectType SNP \
  -o ${file_id}_snp.vcf
"""
}

process filter_indels {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(variants)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_indel.vcf"), emit: vcf
  script:
"""
gatk3 -T SelectVariants \
  -nct ${task.cpus} \
  -R ${fasta} \
  -V ${variants} \
  -selectType INDEL \
  -o ${file_id}_indel.vcf
"""
}

high_confidence_snp_filter = "(QD < 2.0) || (FS > 60.0) || (MQ < 40.0) || (MQRankSum < -12.5) || (ReadPosRankSum < -8.0) || (SOR > 4.0)"

process high_confidence_snp {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(variants)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_snp.vcf"), emit: vcf
  script:
"""
gatk3 -T VariantFiltration \
  -nct ${task.cpus} \
  -R ${fasta} \
  -V ${variants} \
  --filterExpression "${high_confidence_snp_filter}" \
  --filterName "basic_snp_filter" \
  -o ${file_id}_filtered_snp.vcf
"""
}

high_confidence_indel_filter = "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0"

process high_confidence_indel {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(variants)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_indel.vcf"), emit: vcf
  script:
"""
gatk3 -T VariantFiltration \
  -nct ${task.cpus} \
  -R ${fasta} \
  -V ${variants} \
  --filterExpression "${high_confidence_indel_filter}" \
  --filterName "basic_indel_filter" \
  -o ${file_id}_filtered_indel.vcf
"""
}
