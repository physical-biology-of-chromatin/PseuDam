#!/bin/sh
function install_env {
  if [ ! -d "${CONDA_PREFIX}/envs/${1}_${2}" ]; then
    conda create --yes --name ${1}_${2} ${3}=${2}
  fi
}

install_env pigz 2.3.4 pigz
install_env tophat 2.1.1 tophat
install_env hisat2 2.0.0 hisat2
install_env hisat2 2.1.0 hisat2
install_env rsem 1.3.1 rsem
install_env rsem 1.3.0 rsem
install_env samblaster 0.1.24 samblaster
install_env nextflow 0.25.1 nextflow
install_env nextflow 19.01.0 nextflow
install_env nextflow 0.32.0 nextflow
install_env nextflow 0.28.2 nextflow
install_env samtools 1.7 samtools
install_env samtools 1.5 samtools
install_env bowtie2 2.3.2 bowtie2
install_env bowtie2 2.3.4.1 bowtie2
install_env sratools 2.8.2 sra-tools
install_env trimmomatic 0.36 trimmomatic
install_env trimmomatic 0.39 trimmomatic
install_env Python 3.6.1 Python
install_env Python 2.7.13 Python
install_env kallisto 0.44.0 kallisto
install_env kallisto 0.43.1 kallisto
install_env music 1.0.0 music
install_env umitools 0.3.4 umitools
install_env umi_tools 1.0.0 umi_tools
install_env fastp 0.19.7 fastp
install_env gatk 3.8 gatk
install_env cutadapt 1.14 cutadapt
install_env cutadapt 2.1 cutadapt
install_env bioawk 1.0 bioawk
install_env canu 1.7 canu
install_env fastqc 0.11.5 fastqc
install_env bedtools 2.25.0 bedtools
install_env macs2 2.1.2 macs2
install_env bcftools 1.7 bcftools
install_env salmon 0.8.2 salmon
install_env urqt d62c1f8 urqt
install_env multiqc 0.9 multiqc
install_env multiqc 1.7 multiqc
install_env multiqc 1.0 multiqc
install_env cdhit 4.6.8 cdhit
install_env deeptools 3.0.2 deeptools
install_env htseq 0.9.1 htseq
install_env htseq 0.11.2 htseq
install_env python 3.7 python
install_env R 3.5.1 R
install_env R 3.4.3 R
install_env R 3.3.1 R
install_env file handle 0.1.1 file handle
install_env ncdu 1.13 ncdu
install_env picard 2.18.11 picard
install_env sambamba 0.6.7 sambamba

