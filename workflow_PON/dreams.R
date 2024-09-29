#!/usr/bin/env Rscript
reticulate::use_condaenv("dreams_M", required = TRUE)
library(dreams)
library(keras)
args = commandArgs(trailingOnly=TRUE)

ref <- args[1]
pon_bam <- args[2]
pon_data <-args[3]
pon_info <-args[4]

training_data = get_training_data(bam_paths = pon_bam, reference_path = ref, mm_rate_max=0.05)

write.table(training_data$data, pon_data, row.names=FALSE, sep=" ", quote = FALSE)
write.table(training_data$info, pon_info, row.names=FALSE, sep=" ", quote = FALSE)

