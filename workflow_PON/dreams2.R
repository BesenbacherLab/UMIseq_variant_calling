#!/usr/bin/env Rscript
reticulate::use_condaenv("dreams_M", required = TRUE)
library(dreams)
library(keras)
args = commandArgs(trailingOnly=TRUE)
data_files = args[1]
info_files = args[2]

read_data = function(data){
  data <-read.table(data, header=TRUE)
}

combine_training_data = function(data_files, info_files) {
  RDS_data = lapply(data_files, read_data)
  data <- do.call(rbind, RDS_data)
  RDS_info = lapply(info_files, read_data)
  info <- do.call(rbind, RDS_info)
  combined_beta = sum(info$total_matches) / (sum(info$total_coverage) - sum(info$total_mismatches))
  info = data.frame(total_matches = sum(info$total_matches),
                    total_coverage = sum(info$total_coverage),
                    total_mismatches = sum(info$total_mismatches), 
                    beta = combined_beta)
  return(list(data = data, info = info))
}

output_data = combine_training_data(data_files, info_files)


model = train_dreams_model(
    output_data$data,
    layers=c(16,8),
    model_features=c('ref', 'strand', 'first_in_pair', 'read_index', 'trinucleotide_ctx', 'fragment_size', 'umi_count', 'n_other_errors', 'local_GC', 'umi_errors', 'n_deletions_in_read'),
    lr=0.01,
    batch_size=64000, 
    epochs=100, 
    model_file_path=args[3], 
    log_file_path=args[4],
    validation_split = 0.2
)
