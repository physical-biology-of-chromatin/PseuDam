#!/bin/sh

# generate training set
ll *.gz | sed 's/.gz//g' | awk '{system("gunzip -c "$9".gz > ~/data/JU28_59vs17_SNP/data/samples/"$9)}'
cd ~/data/JU28_59vs17_SNP/data/samples/
~/scripts/fastq_sampler/fastq_sampler.py -i 100000 -f MR_350_clean_1.fastq -g MR_350_clean_2.fastq
~/scripts/fastq_sampler/fastq_sampler.py -i 100000 -f MR_550_clean_1.fastq -g MR_550_clean_2.fastq
~/scripts/fastq_sampler/fastq_sampler.py -i 100000 -f NG-10944_JU2859_bis_lib169352_5217_1_1.fastq -g NG-10944_JU2859_bis_lib169352_5217_1_2.fastq
ll s_*.fastq | awk '{system("gzip < "$9" > "$9".gz")}'
cd ~/projects/JU28_59vs17_SNP/

# training set analysis

mkdir tests
cd tests
../nextflow ../src/SNP_calling.nf -c ../src/SNP_calling.config -profile docker --fasta "../data/fasta/DBG2OLC-output2.fasta" --fastq "../data/samples/*_{1,2}.fastq.gz" -resume -w ~/data/work_s/ --tumor "[\"s_NG-10944_JU2859_bis_lib169352_5217_1\"]" --normal "[\"s_MR_550_clean\", \"s_MR_350_clean\"]"
~/scripts/sms.sh "SNP done"

# real set analysis

./nextflow src/SNP_calling.nf -c src/SNP_calling.config -profile docker --fasta "data/fasta/DBG2OLC-output2.fasta" --fastq "data/fastq/*_{1,2}.fastq.gz" -resume -w ~/data/work/ --tumor "[\"NG-10944_JU2859_bis_lib169352_5217_1\"]" --normal "[\"MR_550_clean\", \"MR_350_clean\"]"
~/scripts/sms.sh "SNP done"

./nextflow src/SNP_calling.nf -c src/SNP_calling.config -profile docker --fasta "data/fasta/DBG2OLC-output2.fasta" --fastq "data/fastq/*_{1,2}.fastq.gz" --sam "results/mapping/sam/*.sam" -resume -w ~/data/work/

./nextflow src/SNP_calling.nf -c src/SNP_calling.config -profile docker --fasta "data/fasta/DBG2OLC-output2.fasta" --fastq "data/fastq/*_{1,2}.fastq.gz" --sam "results/mapping/sam/*.sam" -resume -w ~/data/work/ --tumor "[\"NG-10944_JU2859_bis_lib169352_5217_1\"]" --normal "[\"MR_550_clean\", \"MR_350_clean\"]"
~/scripts/sms.sh "SNP done"


