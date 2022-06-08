std_length=150
mean_length=10



files_bed=("/datas/nathan/vscode_nextflow/nextflow-nathan/data/sites/binned_10_100.bed" \
"/datas/nathan/vscode_nextflow/nextflow-nathan/data/sites/binned_10_500.bed" \
"/datas/nathan/vscode_nextflow/nextflow-nathan/data/sites/binned_10_1000.bed" \
"/datas/nathan/vscode_nextflow/nextflow-nathan/data/sites/binned_10_2000.bed" \
"/datas/nathan/vscode_nextflow/nextflow-nathan/data/sites/binned_50_100.bed" \
"/datas/nathan/vscode_nextflow/nextflow-nathan/data/sites/binned_50_500.bed" \
"/datas/nathan/vscode_nextflow/nextflow-nathan/data/sites/binned_50_1000.bed" \
"/datas/nathan/vscode_nextflow/nextflow-nathan/data/sites/binned_50_2000.bed" \
"/datas/nathan/vscode_nextflow/nextflow-nathan/data/sites/binned_100_500.bed" \
"/datas/nathan/vscode_nextflow/nextflow-nathan/data/sites/binned_100_1000.bed" \
"/datas/nathan/vscode_nextflow/nextflow-nathan/data/sites/binned_100_2000.bed" \
"/datas/nathan/vscode_nextflow/nextflow-nathan/data/sites/binned_200_500.bed" \
"/datas/nathan/vscode_nextflow/nextflow-nathan/data/sites/binned_200_1000.bed" \
"/datas/nathan/vscode_nextflow/nextflow-nathan/data/sites/binned_200_2000.bed")



bin_names=("10_100" \
"10_500" \
"10_1000" \
"10_2000" \
"50_100" \
"50_500" \
"50_1000" \
"50_2000" \
"100_500" \
"100_1000" \
"100_2000" \
"200_500" \
"200_1000" \
"200_2000")


for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13
do
./nextflow src/Dam_ID_analysis_v2_bins.nf \
-profile docker \
--mean_length ${mean_length} \
--std_length ${std_length} \
--bed ${files_bed[$i]} \
--bin_dir ${bin_names[$i]} \
-resume
done

