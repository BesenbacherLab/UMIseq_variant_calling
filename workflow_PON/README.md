# UMIseq_variant_calling
## PON processing
**About workflow_PON**

All files under the folder <workflow_PON> are used to process PON files for variant calling by mutect2, shearwater and DREAMS-vc.

**About param.json**

Inside param.json, the working directory, conda environments, gwf workflow, addtional scripts for running, the reference files, the patterns for FASRQ files, and the parameters for UMI grouping and UMI consensus have been stated. 

**Specfic steps** 

***STEP1:*** 

Run the gwf workflow with the templates.py and workflow.py under workflow_PON, to complete the exactly same UMI processing procedures with the plasma samples from the patients. Detailed steps include merging FASTQ files, mapping to HG38, sorting and indexing the BAM file, UMI grouping by umitools, tagging, UMI grouping by fgbio based on the new tags, UMI consensus by fgbio, realignment to HG38, Softclipping the overlapping region, subsetting the target panel, read mate fixation, and generating PILEUP file for each sample.

***STEP2:*** 

The last step included in the gwf workflow is the pon data information required by DREAMS-vc. If you wanna run DREAMS-vc for variant calling, then please run it.

***STEP3:*** 

For DREAMS-vc, an additiona step of model training should be performed on the PON data, and is conducted by dreams2.R under workflow_PON. After changing the directories inside dreams2.R to your own locations, you could run it by `Rscript dreams2.R`.

***STEP4:*** 

For shearwater, it requires to extract PON information as its input, thus, running `Rscript pon_shearwater.R` is necessary if you wanna run shearwater for variant calling.

***STEP5:*** 

For mutect2, if you wanna include PON for filtering, then please run `bash create_pon_vcf.sh` by changing all paths to your own paths.
