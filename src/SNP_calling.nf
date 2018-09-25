params.fastq = "$baseDir/data/*.fastq"
params.fasta = "$baseDir/data/*.fasta"
params.sam = ""
log.info "fastq files : ${params.fastq}"
log.info "fasta files : ${params.fasta}"
def normal_sample = Eval.me(params.normal)
def tumor_sample = Eval.me(params.tumor)
log.info "normal : ${normal_sample}"
log.info "tumor : ${tumor_sample}"

Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .into { fasta_file;
     indel_fasta_file;
     recalibration_fasta_file;
     haplotypecaller_fasta_file
  }
Channel
  .fromFilePairs( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

if (params.sam == "") {
  process adaptor_removal {
    tag "$pair_id"
    publishDir "results/fastq/adaptor_removal/", mode: 'copy'

    input:
    set pair_id, file(reads) from fastq_files

    output:
    set pair_id, "*_cut_R{1,2}.fastq.gz" into fastq_files_cut

    script:
  """
  cutadapt -a AGATCGGAAGAG -g CTCTTCCGATCT -A AGATCGGAAGAG -G CTCTTCCGATCT \
  -o ${pair_id}_cut_R1.fastq.gz -p ${pair_id}_cut_R2.fastq.gz \
  ${reads[0]} ${reads[1]} > ${pair_id}_report.txt
  """
  }

  process trimming {
    tag "${reads}"
    cpus 4
    publishDir "results/fastq/trimming/", mode: 'copy'

    input:
    set pair_id, file(reads) from fastq_files_cut

    output:
    set pair_id, "*_trim_R{1,2}.fastq.gz" into fastq_files_trim

    script:
  """
  UrQt --t 20 --m ${task.cpus} --gz \
  --in ${reads[0]} --inpair ${reads[1]} \
  --out ${pair_id}_trim_R1.fastq.gz --outpair ${pair_id}_trim_R2.fastq.gz \
  > ${pair_id}_trimming_report.txt
  """
  }

  process index_fasta {
    tag "$fasta_id"
    cpus 4
    publishDir "results/mapping/index/", mode: 'copy'

    input:
      set fasta_id, file(fasta) from fasta_file

    output:
      set fasta_id, "${fasta.baseName}.*" into index_files
      file "*_bwa_report.txt" into index_files_report

    script:
  """
  bwa index -p ${fasta.baseName} ${fasta} \
  &> ${fasta.baseName}_bwa_report.txt
  """
  }


  process mapping_fastq {
    tag "$reads"
    cpus 4
    publishDir "results/mapping/sam/", mode: 'copy'

    input:
    set pair_id, file(reads) from fastq_files_trim
    set index_id, file(index) from index_files.collect()

    output:
    file "${pair_id}.sam" into sam_files
    file "${pair_id}_bwa_report.txt" into mapping_repport_files

    script:
  """
  bwa mem -t ${task.cpus} \
  ${index_id} ${reads[0]} ${reads[1]} \
  -o ${pair_id}.sam &> ${pair_id}_bwa_report.txt
  """
  }
} else {
  Channel
    .fromPath( params.sam )
    .ifEmpty { error "Cannot find any sam files matching: ${params.sam}" }
    .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
    .set { sam_files }
}

process dedup_sam {
  tag "$file_id"
  cpus 4

  input:
    set file_id, file(sam) from sam_files

  output:
    set file_id, "*_dedup.sam*" into dedup_sam_files
  script:
"""
samblaster --addMateTags -i ${sam} -o ${file_id}_dedup.sam
"""
}

process sam_to_bam {
  tag "$file_id"
  cpus 4

  input:
    set file_id, file(sam) from dedup_sam_files

  output:
    set file_id, "*.bam" into dedup_bam_files

  script:
"""
sambamba view -t ${task.cpus} -S -f bam -l 0 ${sam} -o ${file_id}.bam
"""
}

process sort_bam {
  tag "$file_id"
  cpus 10

  input:
    set file_id, file(bam) from dedup_bam_files

  output:
    set file_id, "*_sorted.bam" into sorted_bam_files

  script:
"""
sambamba sort -t ${task.cpus} --tmpdir=./tmp -o ${file_id}_sorted.bam ${bam}
"""
}

sorted_bam_files.into {
  sorted_bam_files_norm;
  sorted_bam_files_tumor
}
collect_sorted_bam_file = sorted_bam_files_norm
  .filter{ normal_sample.contains(it[0]) }
  .map { it -> it[1]}
  .collect()
  .map { it -> ["normal_sample", it]}
collect_sorted_bam_file.join(
  sorted_bam_files_tumor
    .filter{ tumor_sample.contains(it[0]) }
    .map { it -> it[1]}
    .collect()
    .map { it -> ["tumor_sample", it]}
)

process merge_bam {
  tag "$file_id"
  cpus 4

  input:
    set file_id, file(bam) from collect_sorted_bam_file

  output:
    set file_id, "*.bam" into merged_bam_files

  script:
"""
sambamba merge -t ${task.cpus} ${file_id}.bam ${bam}
"""
}

process name_bam {
  tag "$file_id"
  cpus 4
  publishDir "results/mapping/bam/", mode: 'copy'

  input:
    set file_id, file(bam) from merged_bam_files

  output:
    set file_id, "*_named.bam" into named_bam_files

  script:
"""
samtools view -H ${bam} > header.sam
echo "@RG\tID:${file_id}\tLB:library1\tPL:illumina\tPU:${file_id}\tSM:${file_id}" \
>> header.sam
cp ${bam} ${file_id}_named.bam
samtools reheader header.sam ${file_id}_named.bam
"""
}

named_bam_files.into{
  index_named_bam_files;
  haplotypecaller_named_bam_files
}

process index_bam {
  tag "$file_id"
  cpus 4
  publishDir "results/mapping/bam/", mode: 'copy'

  input:
    set file_id, file(bam) from index_named_bam_files

  output:
    set file_id, "*.bam*" into indexed_bam_files

  script:
"""
sambamba index -t ${task.cpus} ${bam}
"""
}

haplotypecaller_fasta_file.into{
    haplo_fasta_file;
    index2_fasta_file
    index3_fasta_file
  }

process index2_fasta {
  tag "$genome_id"
  publishDir "results/fasta/", mode: 'copy'

  input:
    set genome_id, file(fasta) from index2_fasta_file

  output:
    set genome_id, "*.dict" into indexed2_fasta_file

  script:
"""
gatk CreateSequenceDictionary -R ${fasta} &> gatk_output.txt
"""
}

process index3_fasta {
  tag "$genome_id"
  publishDir "results/fasta/", mode: 'copy'

  input:
    set genome_id, file(fasta) from index3_fasta_file

  output:
    set genome_id, "*.fai" into indexed3_fasta_file

  script:
"""
samtools faidx ${fasta}
"""
}

process HaplotypeCaller {
  tag "$file_id"
  cpus 4
  publishDir "results/SNP/vcf/", mode: 'copy'

  input:
    set file_id, file(bam) from haplotypecaller_named_bam_files.collect()
    set file_ididx, file(bamidx) from indexed_bam_files.collect()
    set genome_id, file(fasta) from haplo_fasta_file.collect()
    set genome2_idx, file(fasta2idx) from indexed2_fasta_file.collect()
    set genome3_idx, file(fasta3idx) from indexed3_fasta_file.collect()

  output:
    set file_id, "*.vcf" into vcf_files
    set file_id, "*.bam" into realigned_bams_files

  script:
"""
gatk Mutect2 --native-pair-hmm-threads ${task.cpus} -R ${fasta} \
-I ${bam} -tumor ${params.tumor} -normal ${params.normal} \
-O ${file_id}_raw_calls.g.vcf \
-bamout ${file_id}_realigned.bam
"""
}

/*
process filter_SNP {
  tag "$file_id"
  cpus 4
  publishDir "results/SNP/vcf/", mode: 'copy'

  input:

  output:
    set file_id, "*.vcf" into vcf_files_filtered

  script:
"""
gatk --java-options "-Xmx2g" Mutect2 \
-R hg38/Homo_sapiens_assembly38.fasta \
-I tumor.bam \
-I normal.bam \
-tumor HCC1143_tumor \
-normal HCC1143_normal \
-pon resources/chr17_pon.vcf.gz \
--germline-resource resources/chr17_af-only-gnomad_grch38.vcf.gz \
--af-of-alleles-not-in-resource 0.0000025 \
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
-L chr17plus.interval_list \
-O 1_somatic_m2.vcf.gz \
-bamout 2_tumor_normal_m2.bam

gatk Mutect2 \
-R ~/Documents/ref/hg38/Homo_sapiens_assembly38.fasta \
-I HG00190.bam \
-tumor HG00190 \
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
-L chr17plus.interval_list \
-O 3_HG00190.vcf.gz

gatk CreateSomaticPanelOfNormals \
-vcfs 3_HG00190.vcf.gz \
-vcfs 4_NA19771.vcf.gz \
-vcfs 5_HG02759.vcf.gz \
-O 6_threesamplepon.vcf.gz

gatk GetPileupSummaries \
-I tumor.bam \
-V resources/chr17_small_exac_common_3_grch38.vcf.gz \
-O 7_tumor_getpileupsummaries.table

gatk CalculateContamination \
-I 7_tumor_getpileupsummaries.table \
-O 8_tumor_calculatecontamination.table

gatk FilterMutectCalls \
-V somatic_m2.vcf.gz \
--contamination-table tumor_calculatecontamination.table \
-O 9_somatic_oncefiltered.vcf.gz

gatk CollectSequencingArtifactMetrics \
-I tumor.bam \
-O 10_tumor_artifact \
–-FILE_EXTENSION ".txt" \
-R ~/Documents/ref/hg38/Homo_sapiens_assembly38.fasta

gatk FilterByOrientationBias \
-A G/T \
-A C/T \
-V 9_somatic_oncefiltered.vcf.gz \
-P tumor_artifact.pre_adapter_detail_metrics.txt \
-O 11_somatic_twicefiltered.vcf.gz

"""
}

*/
