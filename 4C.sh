for chrom in I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI
do
    hicPlotMatrix --matrix "/home/nathan/Downloads/GSE151553_A364_merged.mcool::resolutions/50" --outFileName 4C_chrom${chrom}.eps --log --clearMaskedBins --region ${chrom}:1-1000000 --region2 III:80000-100000
done