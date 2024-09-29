# UMIseq_variant_calling
## PON processing
**About <workflow_PON>**

All files under the folder <workflow_PON> are used to process PON files and prepare all files in need for variant calling by mutect2, shearwater and DREAMS-vc.

**About param.json**

Inside param.json, the working directory, conda environments, gwf workflow, addtional scripts for running, the reference files, the patterns for FASRQ files, and the parameters for UMI grouping and UMI consensus have been stated. 

**Specfic steps** 

***STEP1:UMI processing*** 

Run the gwf workflow with the templates.py and workflow.py under workflow_PON, to complete the exactly same UMI processing procedures with the plasma samples from the patients. Detailed steps include merging FASTQ files, mapping to hg38, sorting and indexing the BAM file, UMI grouping by umitools, tagging, UMI grouping by fgbio based on the new tags, UMI consensus by fgbio, realignment to hg38, softclipping the overlapping region, subsetting the target panel, read mate fixation, and generating PILEUP file for each sample.

***STEP2: Prepare for variant calling*** 

All steps after PILEUP generation in the gwf workflow are performed to prepare files required by DREAMS-vc, shearwater and Mutect2. If you wanna run these variant callers for variant calling, then please also run these steps.
