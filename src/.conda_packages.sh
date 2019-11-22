source src/.conda_psmn.sh
CONDA_ENVS=src/.conda_envs/
if [ ! -d ${CONDA_ENVS}pigz_2.3.4 ]; then
  conda create --yes --name pigz_2.3.4 pigz=2.3.4
fi
if [ ! -d ${CONDA_ENVS}tophat_2.1.1 ]; then
  conda create --yes --name tophat_2.1.1 tophat=2.1.1
fi
if [ ! -d ${CONDA_ENVS}hisat2_2.0.0 ]; then
  conda create --yes --name hisat2_2.0.0 hisat2=2.0.0 samtools=1.7
fi
if [ ! -d ${CONDA_ENVS}hisat2_2.1.0 ]; then
  conda create --yes --name hisat2_2.1.0 hisat2=2.1.0 samtools=1.7
fi
if [ ! -d ${CONDA_ENVS}rsem_1.3.1 ]; then
  conda create --yes --name rsem_1.3.1 rsem=1.3.1 samtools=1.3
fi
if [ ! -d ${CONDA_ENVS}rsem_1.3.0 ]; then
  conda create --yes --name rsem_1.3.0 rsem=1.3.0 samtools=1.3
fi
if [ ! -d ${CONDA_ENVS}samblaster_0.1.24 ]; then
  conda create --yes --name samblaster_0.1.24 samblaster=0.1.24
fi
if [ ! -d ${CONDA_ENVS}nextflow_0.25.1 ]; then
  conda create --yes --name nextflow_0.25.1 nextflow=0.25.1
fi
if [ ! -d ${CONDA_ENVS}nextflow_19.01.0 ]; then
  conda create --yes --name nextflow_19.01.0 nextflow=19.01.0
fi
if [ ! -d ${CONDA_ENVS}nextflow_0.32.0 ]; then
  conda create --yes --name nextflow_0.32.0 nextflow=0.32.0
fi
if [ ! -d ${CONDA_ENVS}nextflow_0.28.2 ]; then
  conda create --yes --name nextflow_0.28.2 nextflow=0.28.2
fi
if [ ! -d ${CONDA_ENVS}samtools_1.7 ]; then
  conda create --yes --name samtools_1.7 samtools=1.7
fi
if [ ! -d ${CONDA_ENVS}samtools_1.5 ]; then
  conda create --yes --name samtools_1.5 samtools=1.5
fi
if [ ! -d ${CONDA_ENVS}bowtie2_2.3.2 ]; then
  conda create --yes --name bowtie2_2.3.2 bowtie2=2.3.2 samtools=1.7
fi
if [ ! -d ${CONDA_ENVS}bowtie2_2.3.4.1 ]; then
  conda create --yes --name bowtie2_2.3.4.1 bowtie2=2.3.4.1 samtools=1.7 #&& \
fi
if [ ! -d ${CONDA_ENVS}sra-tools_2.8.2 ]; then
  conda create --yes --name sra-tools_2.8.2 sra-tools=2.8.2
fi
if [ ! -d ${CONDA_ENVS}trimmomatic_0.36 ]; then
  conda create --yes --name trimmomatic_0.36 trimmomatic=0.36
fi
if [ ! -d ${CONDA_ENVS}trimmomatic_0.39 ]; then
  conda create --yes --name trimmomatic_0.39 trimmomatic=0.39
fi
if [ ! -d ${CONDA_ENVS}Python_3.6.1 ]; then
  conda create --yes --name Python_3.6.1 Python=3.6.1
fi
if [ ! -d ${CONDA_ENVS}Python_2.7.13 ]; then
  conda create --yes --name Python_2.7.13 Python=2.7.13
fi
if [ ! -d ${CONDA_ENVS}kallisto_0.44.0 ]; then
  conda create --yes --name kallisto_0.44.0 kallisto=0.44.0
fi
if [ ! -d ${CONDA_ENVS}kallisto_0.43.1 ]; then
  conda create --yes --name kallisto_0.43.1 kallisto=0.43.1
fi
if [ ! -d ${CONDA_ENVS}music_1.0.0 ]; then
  conda create --yes --name music_1.0.0 music=1.0.0
fi
if [ ! -d ${CONDA_ENVS}umitools_0.3.4 ]; then
  conda create --yes --name umitools_0.3.4 umitools=0.3.4
fi
if [ ! -d ${CONDA_ENVS}fastp_0.19.7 ]; then
  conda create --yes --name fastp_0.19.7 fastp=0.19.7
fi
if [ ! -d ${CONDA_ENVS}gatk_3.8 ]; then
  conda create --yes --name gatk_3.8 gatk=3.8
fi
if [ ! -d ${CONDA_ENVS}cutadapt_1.14 ]; then
  conda create --yes --name cutadapt_1.14 cutadapt=1.14
fi
if [ ! -d ${CONDA_ENVS}bioawk_1.0 ]; then
  conda create --yes --name bioawk_1.0 bioawk=1.0
fi
if [ ! -d ${CONDA_ENVS}canu_1.7 ]; then
  conda create --yes --name canu_1.7 canu=1.7
fi
if [ ! -d ${CONDA_ENVS}fastqc_0.11.5 ]; then
  conda create --yes --name fastqc_0.11.5 fastqc=0.11.5
fi
if [ ! -d ${CONDA_ENVS}bedtools_2.25.0 ]; then
  conda create --yes --name bedtools_2.25.0 bedtools=2.25.0
fi
if [ ! -d ${CONDA_ENVS}macs2_2.1.2 ]; then
  conda create --yes --name macs2_2.1.2 macs2=2.1.2
fi
if [ ! -d ${CONDA_ENVS}bcftools_1.7 ]; then
  conda create --yes --name bcftools_1.7 bcftools=1.7
fi
if [ ! -d ${CONDA_ENVS}salmon_0.8.2 ]; then
  conda create --yes --name salmon_0.8.2 salmon=0.8.2
fi
if [ ! -d ${CONDA_ENVS}urqt_d62c1f8 ]; then
  conda create --yes --name urqt_d62c1f8 urqt=d62c1f8
fi
if [ ! -d ${CONDA_ENVS}multiqc_0.9 ]; then
  conda create --yes --name multiqc_0.9 multiqc=0.9
fi
if [ ! -d ${CONDA_ENVS}multiqc_1.7 ]; then
  conda create --yes --name multiqc_1.7 multiqc=1.7
fi
if [ ! -d ${CONDA_ENVS}multiqc_1.0 ]; then
  conda create --yes --name multiqc_1.0 multiqc=1.0
fi
if [ ! -d ${CONDA_ENVS}cdhit_4.6.8 ]; then
  conda create --yes --name cdhit_4.6.8 cdhit=4.6.8
fi
if [ ! -d ${CONDA_ENVS}deeptools_3.0.2 ]; then
  conda create --yes --name deeptools_3.0.2 deeptools=3.0.2
fi
if [ ! -d ${CONDA_ENVS}htseq_0.9.1 ]; then
  conda create --yes --name htseq_0.9.1 htseq=0.9.1
fi
if [ ! -d ${CONDA_ENVS}htseq_0.11.2 ]; then
  conda create --yes --name htseq_0.11.2 htseq=0.11.2
fi
if [ ! -d ${CONDA_ENVS}R_3.4.3 ]; then
  conda create --yes --name R_3.4.3 R=3.4.3
fi
if [ ! -d ${CONDA_ENVS}R_3.3.1 ]; then
  conda create --yes --name R_3.3.1 R=3.3.1
fi
if [ ! -d ${CONDA_ENVS}file_handle_0.1.1 ]; then
  conda create --yes --name file_handle_0.1.1 file_handle=0.1.1
fi
if [ ! -d ${CONDA_ENVS}ncdu_1.13 ]; then
  conda create --yes --name ncdu_1.13 ncdu=1.13
fi
if [ ! -d ${CONDA_ENVS}picard_2.18.11 ]; then
  conda create --yes --name picard_2.18.11 picard=2.18.11
fi
if [ ! -d ${CONDA_ENVS}sambamba_0.6.7 ]; then
  conda create --yes --name sambamba_0.6.7 sambamba=0.6.7
fi
if [ ! -d ${CONDA_ENVS}star_2.7.3a ]; then
  conda create --yes --name star_2.7.3a star=2.7.3a
fi
if [ ! -d ${CONDA_ENVS}liftover_357 ]; then
  conda create --yes --name liftover_357 ucsc-liftover==357
fi
if [ ! -d ${CONDA_ENVS}axtchain_377 ]; then
  conda create --yes --name axtchain_377 ucsc-axtchain==377
fi
if [ ! -d ${CONDA_ENVS}ribotish_0.2.4 ]; then
  conda create --name ribotish_0.2.4 ribotish=0.2.4
fi
