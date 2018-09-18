#!/bin/sh
./nextflow src/SNP_calling.nf -c src/SNP_calling.config -profile docker --fasta "data/fasta/DBG2OLC-output2.fasta" --fastq "data/fastq/*_{1,2}.fastq.gz" -resume -w ~/data/work/

