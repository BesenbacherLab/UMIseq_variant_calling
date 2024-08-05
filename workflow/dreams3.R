#!/usr/bin/env Rscript
reticulate::use_condaenv("dreams_M", required = TRUE)
library(dreams)
library(keras)
args = commandArgs(trailingOnly=TRUE)
ref <- args[1]
mutations <-read.table(args[2], header=TRUE)
bam <- args[3]
model <- load_model_hdf5(args[4])
beta <- read.table(args[5], header=TRUE)$beta
bam_calls <-args[6]

#####
mutation_calls <-dreams_vc(
  mutations_df = mutations, 
  bam_file_path = bam, 
  reference_path = ref, 
  model = model)

write.table(mutation_calls, bam_calls, row.names=FALSE, sep=" ", quote = FALSE)
