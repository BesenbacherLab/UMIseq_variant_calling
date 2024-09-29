# UMIseq_variant_calling
The codes are used for processing deep-targeted UMI-seq data. 

The workflows are written in gwf (<https://gwf.app>).

It includes the steps of both UMI processing and variant calling options of 

shearwater(both AND & OR algorithm) (<https://github.com/im3sanger/deepSNV>), 

Mutect2 (<https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2>), 

Varscan2 (<https://varscan.sourceforge.net>), 

DREAMS-vc (<https://github.com/JakobSkouPedersenLab/dreams>). 

We used it for benchmarking different variant callers and UMI processing strategies for detecting low-frequency mutations from cfDNA colorectal cancer patients. For detailed theories behind each method, and the specific step for modification based on your own use, please also check the related github respectively.

## Directory ***workflow*** 

Directory ***workflow*** includes all steps for processing cancer patient plasma sample, from UMI processing to variant calling. Please check the README under it for more information. 

For running GWF workflow, please locate at the working directory which includes both the `templates.py` and `workflow.py`. Specifically, the steps under ***workflow***, can be performed by

```
cd workflow
gwf -b slurm run
``` 

Remember to change the **account, running time, core and memory** for the functions inside templates.py for your own use.

## Directory ***workflow_PON*** 

Directory ***workflow_PON*** includes all steps for processing PON which is required by variant calling of mutect2, shearwater and DREAMS-vc. Therefore, before performing variant calling for cancer plasma sample, please check the README under this directory and follow all steps to prepare all files in need. 

For running GWF workflow, please locate at the working directory which includes both the `templates.py` and `workflow.py`. 

Specifically, the steps under ***workflow_PON***, can be performed by  

```
cd workflow_pon
gwf -b slurm run
``` 

Remember to change the **account, running time, core and memory** for the functions inside templates.py for your own use.

## Conda environments

As shown by the json file under workflow and workflow_PON, we utilized three different conda environments to run the workflow. 

1. For conda environment named 'umiseq', please run 

```
conda create -n umiseq -c bioconda -c gwforg python=3.7 gwf samtools bwa picard umi_tools fgbio gatk4 pysamstats seqtk
```

please check `environment_umiseq.yml` for all dependencies and related packages. 

2. For conda environment named 'shearwater', please run

```
conda create -n shearwater
conda activate shearwater
conda install conda-forge::r-tidyverse
```

To install package 'deepSNV', start R (version "4.4") and enter:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("deepSNV")```

please check `environment_sw.yml` for all dependencies and related packages.

3. For conda environment named 'dreams_M', please run

```
conda create -n dreams_M
conda activate dreams_M
```

Under R, install DREAMS package by

```
install.packages("devtools")
devtools::install_github("JakobSkouPedersenLab/dreams")
```

please check `environment_dreams.yml` for all dependencies and related packages. 





