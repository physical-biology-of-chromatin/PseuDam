#!/bin/bash
# TODO remove when process is done

input="hello"
output=$(python3 reads_infos.py --fastq /datas/nathan/vscode_nextflow/nextflow-nathan/data/reads/Dam_ID/c_elegans/SRR13429076.fastq)
echo $output
echo $input