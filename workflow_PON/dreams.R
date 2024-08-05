"
\033[1;4;32m
Dreams read pon data
\033[0m

Usage:
  dreams.R [-h] --ref=REF --pondata=PONDATA --poninfo=PONINFO <ponbam>
  
Options:
-h --help                 Show this screen.
--ref=REF                 Reference fasta.
--pondata=PONDATA           DATA output after reading pon_bam data.
--poninfo=PONINFO           INFO output after reading pon_bam data.

<ponbam>         Pon sorted_targeted_bam file.


Examples:
dreams.R     --ref='ref.fasta' --pondata='pon_data.csv' --poninfo='pon_info.csv' 'pon_bams/Donor284_sorted_target.bam'
Requires:
docopt, dreams

" -> doc
suppressWarnings(suppressMessages(require(docopt)))
suppressWarnings(suppressMessages(require(dreams)))


opt <- docopt(doc)

cat("[dreams.R] Running with arguments:\n")
print(opt)
cat("\n")

ref <- opt$ref
pon_bam <- opt$ponbam
pon_data <-opt$pondata
pon_info <-opt$poninfo

training_data = get_training_data(bam_paths = pon_bam, reference_path = ref, mm_rate_max=0.05)

write.table(training_data$data, pon_data, row.names=FALSE, sep=" ", quote = FALSE)
write.table(training_data$info, pon_info, row.names=FALSE, sep=" ", quote = FALSE)

