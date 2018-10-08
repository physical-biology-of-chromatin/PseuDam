#!/usr/bin/Rscript
library("tidyverse")

args <- commandArgs(trailingOnly = TRUE)
snp_a <- read_delim(args[1], delim = "\t") %>%
  mutate(cords = paste0(CHROM, POS))
snp_b <- read_delim(args[2], delim = "\t") %>%
  mutate(cords = paste0(CHROM, POS))


only_b <- snp_b %>%
  select(cords) %>%
  setdiff(snp_a %>% select(cords)) %>%
  pull(cords)

snp_b %>%
  filter(cords %in% only_b) %>%
  write_csv(paste0("only_", args[2]))

