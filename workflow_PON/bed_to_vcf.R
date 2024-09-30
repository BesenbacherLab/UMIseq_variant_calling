#!/usr/bin/env Rscript
##For calling at all positions in a BED-file you will have to create a vcf-file containing all the positions. 
library('tidyverse')
# Read panel from bed file
args = commandArgs(trailingOnly=TRUE)
panel_bed = read.table(arg[1], sep = "\t", col.names = c("CHROM", "start", "stop", "region")) %>% 
  select(-region) %>% 
  mutate(start = start, stop = stop)
 
# Extend to get df per position
mutations  <-
  panel_bed %>%
  rowwise() %>%
  mutate(POS = list(start:stop)) %>%
  unnest(cols = c("POS")) %>%
  select(-start, -stop)
 
# Initiate vcf file
vcf_out <-file(arg[2])
writeLines(c("##fileformat=VCFv4.1",
             "##fileDate=20220615",
             "##source=myImputationProgramV3.1",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"),
           vcf_out)
 
 
# Prepare mutation candidates for vcf format
pileup = read.table(arg[3], sep="\t", header=TRUE)

vcf = mutations %>% mutate(ID = ".", QUAL = 40, FILTER = ".", INFO = ".")
 
vcf_1 = vcf %>% filter(!CHROM %in% c("chrX", "chrY")) %>% mutate(CHROM_idx = as.integer(str_remove(CHROM, "chr"))) %>% arrange(CHROM_idx, POS) %>% select(CHROM, POS, ID, QUAL, FILTER, INFO)
 
vcf_2 = vcf %>% filter(CHROM %in% c("chrX", "chrY")) %>% mutate(CHROM_idx = str_remove(CHROM, "chr")) %>% arrange(CHROM_idx, POS) %>% select(CHROM, POS, ID, QUAL, FILTER, INFO)
 
vcf_out = rbind(vcf_1, vcf_2) %>% distinct()
vcf_ref = left_join(vcf_out, pileup[,1:3], by=c("CHROM"="chrom", "POS"="pos")) %>% relocate(REF = ref, .after = ID)
vcf_alt <-vcf_ref %>%
  rowwise() %>%
  mutate(ALT = list(c("A", "T", "C", "G"))) %>%
  unnest(cols = c("ALT")) %>%
  relocate(ALT, .after = REF) %>%
  filter(REF!=ALT)

# Write to vcf file
#for mutect2
write.table(vcf_alt, arg[2], append = T, sep = "\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
#for dreams-vc
write.table(vcf_alt, arg[4], sep = "\t", quote=FALSE, col.names=TRUE, row.names=FALSE)