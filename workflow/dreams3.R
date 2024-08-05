#"
#\033[1;4;32m
#Dreams call variants
#\033[0m
#
#Usage:
#  dreams3.R [-h] --mutations=MUTATIONS --ref=REF --model=MODEL --allponinfo=ALLPONINFO --calls=CALLS --batchsize=BATCHSIZE --ncores=NCORES --logfile=LOGFILE <bam>
#  
#Options:
#-h --help                 Show this screen.
#--mutations=MUTATIONS     Mutations forced to call.
#--ref=REF                 Reference fasta.
#--model=MODEL             Model after training.
#--allponinfo=ALLPONINFO   To get Beta.
#--calls=CALLS             Output variant calls.
#--batchsize=BATCHSIZE     Number of position to run at the same time.
#--ncores=NCORES           Number of cores
#--logfile=LOGFILE         Log file
#
#<bam>         Sorted_targeted_bam file to call variants.
#
#
#Examples:
#dreams3.R     --ref='ref.fasta' --allponinfo='all_pon_info.csv' --mutations='mutations.df' --model='model.hdf5' --calls='fq_data2/B05464/B05464_dreams_variants.df' --batchsize=2 --ncores=20 --logfile='fq_data2/B05464/B05464_dreams.log' 'fq_data2/B05464/B05464_sorted_target_consensus.bam'
#Requires:
#docopt, dreams
#
#" -> doc
#suppressWarnings(suppressMessages(require(docopt)))
#suppressWarnings(suppressMessages(require(dreams)))


#opt <- docopt(doc)
#
#cat("[dreams3.R] Running with arguments:\n")
#print(opt)
#cat("\n")

#ref <- opt$ref
#mutations <-read.table(opt$mutations, header=TRUE)
#bam <- opt$bam
#model <- keras::load_model_hdf5(opt$model)
#beta <- read.table(opt$allponinfo, header=TRUE)$beta
#bam_calls <-opt$calls
#batch_size <-as.numeric(opt$batchsize)
#ncores <-as.numeric(opt$ncores)
#log_file <-opt$logfile

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
#batch_size <-2
#ncores <-15
#log_file <-args[7]

#####
mutation_calls <-dreams_vc( #_parallel
  mutations_df = mutations, 
  bam_file_path = bam, 
  reference_path = ref, 
  model = model)#, 
  #beta = beta, 
  #batch_size = batch_size, 
  #ncores = ncores,
  #log_file=log_file) 

write.table(mutation_calls, bam_calls, row.names=FALSE, sep=" ", quote = FALSE)
