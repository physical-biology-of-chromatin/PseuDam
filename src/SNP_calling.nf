params.fastq = "$baseDir/data/*.fastq"
params.fasta = "$baseDir/data/*.fasta"
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

process adaptor_removal {
  tag "$pair_id"
  publishDir "results/fastq/adaptor_removal/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files

  output:
  set pair_id, file("*.fastq.gz") into fastq_files_cut
  file "*_cutadapt_report.txt" into cut_files_report

  script:
"""
cutadapt -a AGATCGGAAGAG -g CTCTTCCGATCT -A AGATCGGAAGAG -G CTCTTCCGATCT \
-o ${pair_id}_cut_R1.fastq.gz -p ${pair_id}_cut_R2.fastq.gz \
${reads[0]} ${reads[1]} > ${pair_id}_cutadapt_report.txt
"""
}

process trimming {
  tag "${pair_id}"
  cpus 4
  publishDir "results/fastq/trimming/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files_cut

  output:
  set pair_id, file("*.fastq.gz") into fastq_files_trim
  file "*_trimming_report.txt" into trimming_files_report

  script:
"""
UrQt --t 20 --m ${task.cpus} --gz \
--in ${reads[0]} --inpair ${reads[1]} \
--out ${pair_id}_trim_R1.fastq.gz --outpair ${pair_id}_trim_R2.fastq.gz \
> ${pair_id}_trimming_report.txt
"""
}

process index_fasta {
  tag "$file_id"
  cpus 4
  publishDir "results/mapping/index/", mode: 'copy'

  input:
    set file_id, file(fasta) from fasta_file

  output:
    file "*.index*" into index_files
    file "*_report.txt" into indexing_report

  script:
"""
bowtie2-build --threads ${task.cpus} ${fasta} ${file_id}.index &> ${file_id}_bowtie2_report.txt

if grep -q "Error" ${file_id}_bowtie2_report.txt; then
  exit 1
fi
"""
}

process mapping_fastq {
  tag "$pair_id"
  cpus 4
  publishDir "results/mapping/bams/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files_trim
  file index from index_files.collect()

  output:
  set pair_id, "*.bam" into bam_files
  file "*_report.txt" into mapping_report

  script:
  index_id = index[0]
  for (index_file in index) {
    if (index_file =~ /.*\.1\.bt2/ && !(index_file =~ /.*\.rev\.1\.bt2/)) {
        index_id = ( index_file =~ /(.*)\.1\.bt2/)[0][1]
    }
  }
"""
bowtie2 --very-sensitive -p ${task.cpus} -x ${index_id} \
-1 ${reads[0]} -2 ${reads[1]} 2> \
${pair_id}_bowtie2_report.txt | \
samtools view -Sb - > ${pair_id}.bam

if grep -q "Error" ${pair_id}_bowtie2_report.txt; then
  exit 1
fi
"""
}

process sort_bam {
  tag "$file_id"
  cpus 4

  input:
    set file_id, file(bam) from bam_files

  output:
    set file_id, "*_sorted.bam" into sorted_bam_files

  script:
"""
sambamba sort -t ${task.cpus} --tmpdir=./tmp -o ${file_id}_sorted.bam ${bam}
"""
}

sorted_bam_files.into {
  sorted_bam_file_norm;
  sorted_bam_file_tumor
}

collect_sorted_bam_file_norm = sorted_bam_file_norm
  .filter{ normal_sample.contains(it[0]) }
  .map { it -> it[1]}
  .collect()
  .map { it -> ["normal_sample", it]}

collect_sorted_bam_file_tumor = sorted_bam_file_tumor
  .filter{ tumor_sample.contains(it[0]) }
  .map { it -> it[1]}
  .collect()
  .map { it -> ["tumor_sample", it]}

collect_sorted_bam_file = Channel.create()
  .mix(collect_sorted_bam_file_norm, collect_sorted_bam_file_tumor)

process merge_bam {
  tag "$file_id"
  cpus 4

  input:
    set file_id, file(bam) from collect_sorted_bam_file

  output:
    set file_id, "*.bam" into merged_bam_files

  script:
"""
if ((\$(ls -l *.bam | wc -l) > 1)); then
sambamba merge -t ${task.cpus} ${file_id}.bam ${bam}
else
cp ${bam} ${file_id}.bam
fi
"""
}

merged_bam_files.into{
  index_merged_bam_files;
  haplo_bam_files_norm;
  haplo_bam_files_tumor
}

process index_bam {
  tag "$file_id"
  cpus 4
  publishDir "results/mapping/bam/", mode: 'copy'

  input:
    set file_id, file(bam) from index_merged_bam_files

  output:
    set file_id, "*.bam.bai" into index_bam_files

  script:
"""
sambamba index -t ${task.cpus} ${bam}
"""
}

index_bam_files.into{
  named_index_bam_files;
  indexed_bam_files
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

haplotypecaller_bam_files_norm = haplo_bam_files_norm
  .filter{ "normal_sample" == it[0] }
haplotypecaller_bam_files_tumor = haplo_bam_files_tumor
   .filter{ "tumor_sample" == it[0] }

indexed_bam_files.into {
  index_bam_files_norm;
  index_bam_files_tumor
}
indexed_bam_files_norm = index_bam_files_norm
  .filter{ "normal_sample" == it[0] }
indexed_bam_files_tumor = index_bam_files_tumor
   .filter{ "tumor_sample" == it[0] }

process HaplotypeCaller {
  tag "$file_id"
  cpus 4
  publishDir "results/SNP/vcf/", mode: 'copy'

  input:
    set file_id_norm, file(bam_norm) from haplotypecaller_bam_files_norm.collect()
    set file_ididx_norm, file(bamidx_norm) from indexed_bam_files_norm.collect()
    set file_id_tumor, file(bam_tumor) from haplotypecaller_bam_files_tumor.collect()
    set file_ididx_tumor, file(bamidx_tumor) from indexed_bam_files_tumor.collect()
    set genome_id, file(fasta) from haplo_fasta_file.collect()
    set genome2_idx, file(fasta2idx) from indexed2_fasta_file.collect()
    set genome3_idx, file(fasta3idx) from indexed3_fasta_file.collect()

  output:
    set file_id, "*.vcf" into vcf_files
    set file_id, "*.bam" into realigned_bams_files

  script:
"""
gatk Mutect2 --native-pair-hmm-threads ${task.cpus} -R ${fasta} \
-I ${bam_tumor} -tumor ${file_id_tumor} \
-I ${bam_norm} -normal ${file_id_norm} \
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
