#!/usr/bin/Rscript
library("tidyverse")
library("seqinr")

args <- commandArgs(trailingOnly = TRUE)
snp_a <- read_delim(args[1], delim = "\t") %>%
  mutate(cords = paste0(CHROM, POS))
snp_b <- read_delim(args[2], delim = "\t") %>%
  mutate(cords = paste0(CHROM, POS))


only_b <- snp_b %>%
  select(cords) %>%
  setdiff(snp_a %>% select(cords)) %>%
  pull(cords)
snp <- snp_b %>%
  filter(cords %in% only_b)

fastafile <- read.fasta(file = args[3],
                        as.string = TRUE)

snp$seq_list <- snp %>%
  apply(1, FUN = function(x, fastafile, POS, CHROM){
      begin <- as.integer(x[ POS ]) - 10
      end <- as.integer(x[ POS ]) + 10
      chrom <- x[ CHROM ]
      seq_restric <- fastafile[[ chrom ]] %>%
        c2s() %>%
        substr(begin, end) %>%
        s2c()
      seq_restric[11] <- toupper(seq_restric[11])
      return(seq_restric %>% c2s())
    },
    fastafile = fastafile,
    POS = which(colnames(snp) %in% "POS"),
    CHROM = which(colnames(snp) %in% "CHROM")
  )

snp %>%
  arrange(desc(tumor_sample.AF)) %>%
  write_csv(paste0("only_", args[2]))

