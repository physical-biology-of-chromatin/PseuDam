#!/bin/sh
./nextflow src/SNP_calling.nf -c src/SNP_calling.config -profile docker --fasta "data/fasta/DBG2OLC-output2.fasta" --fastq "data/fastq/*_{1,2}.fastq.gz" -resume -w ~/data/work/

./nextflow src/SNP_calling.nf -c src/SNP_calling.config -profile docker --fasta "data/fasta/DBG2OLC-output2.fasta" --fastq "data/fastq/*_{1,2}.fastq.gz" --sam "results/mapping/sam/*.sam" -resume -w ~/data/work/

./nextflow src/SNP_calling.nf -c src/SNP_calling.config -profile docker --fasta "data/fasta/DBG2OLC-output2.fasta" --fastq "data/fastq/*_{1,2}.fastq.gz" --sam "results/mapping/sam/*.sam" -resume -w ~/data/work/ --tumor "NG-10944_JU2859_bis_lib169352_5217_1" --normal "MR_550_clean"
