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

## Directory ***workflow_PON*** 

Directory ***workflow_PON*** includes all steps for processing PON which is required by variant calling of mutect2, shearwater and DREAMS-vc. Therefore, before performing variant calling for cancer plasma sample, please check the README under this directory and follow all steps to prepare all files in need. 

## GWF workflow

For running GWF workflow, please locate at the working directory which includes both the `templates.py` and `workflow.py`. It is easy to run by `gwf -b slurm run`. Also, remember to change the **account, running time, core and memory** for the functions inside templates.py.

## Conda environments

As shown by the json file under workflow and workflow_PON, we utilized three different conda environments to run the workflow. 

1. For conda environment named 'umiseq', please check `environment_umiseq.yml` for all dependencies and related packages. 

2. For conda environment named 'shearwater', please check `environment_sw.yml` for all dependencies and related packages.

3. For conda environment named 'dreams_M', please check `environment_dreams.yml` for all dependencies and related packages. 





