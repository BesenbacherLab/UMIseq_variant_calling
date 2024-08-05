import os

#inputs = "/home/yixinlin/ctdna_var_calling/PrimaryData"
#outputs = "/home/yixinlin/ctdna_var_calling/Workspaces/yixinlin/fq_data2"
#Change to your directory

samples = os.listdir(inputs)
for sample in samples:
    if os.path.isfile(os.path.join(inputs, sample)):
       samples.remove(sample)
for sample in samples:
    os.mkdir(os.path.join(outputs, sample))

def create_fq_list(sample):
    """
    for each sample, create a fq_list for fastq files, where each line
    is a lane, each column is a read.
    (e.g. a_R1_001.fastq.gz a_R2_001_fastq.gz a_UMI_001_fastq.gz\n b_R1...)
    """

    fq_list = os.listdir(os.path.join(inputs, sample))
    true_fq_list = []
    for fq in fq_list:
        if fq.endswith(("R1_001.fastq.gz", "R2_001.fastq.gz", "UMI_001.fastq.gz")):
            true_fq_list.append(os.path.join(inputs, sample, fq))
    textfile = open(os.path.join(outputs, f"{sample}/{sample}_fq.lst"), "w")
    for i, n in enumerate(sorted(true_fq_list)):
        textfile.write(n + " ")
        if i%3 ==2:
            textfile.write("\n")
    textfile.close()
    return textfile


for sample in samples:
    create_fq_list(sample)


