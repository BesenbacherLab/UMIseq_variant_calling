# UMIseq_variant_calling
## Patient plasma sample processing
**About workflow**

All files under the folder <workflow> are used to process patient plasma sample files. It includes the UMI processing step and variant calling by VarScan2, Mutect2, shearwater and DREAMS-vc.

**About param.json**

Inside param.json, the working directory, conda environments, gwf workflow, addtional scripts for running, the reference files, the patterns for FASRQ files, and the parameters for UMI grouping and UMI consensus have been stated. Please check whether all the reference files exist before you run the workflow, and change the file paths and working directory for your own use. 

**Specfic steps** 

***STEP1: UMI processing*** 

Run the gwf workflow with the templates.py and workflow.py under workflow to process each patient plasma sample file. Detailed steps include merging FASTQ files, mapping to HG38, sorting and indexing the BAM file, UMI grouping by umitools, tagging, UMI grouping by fgbio based on the new tags, UMI consensus by fgbio, realignment to HG38, Softclipping the overlapping region, subsetting the target panel, read mate fixation, and generating PILEUP file for each sample.

In addition to the UMI processing steps above, variant calling steps are also included in the code. For filtering non-cancer mutations and SNPs, PILEUP and MPILEUP are also generated for buffycoat(or saying PBMC) sample.

***STEP2: Run shearwater***
Before running the gwf workflow step of shearwtaer variant calling, please check if `pon_sw.RDS` exists after running all steps shown under <workflow_PON>. Then run variant calling by shearwater following the gwf workflow. 

***STEP3: Run Mutect2***

For running mutect2 under tumor-normal mode and output all mutations along the panel, additional step is needed before variant calling. Please run `Rscript bed_vcf_mutect2.R` to output a dataframe (named `mutect2_allbedpos.vcf`) including all positions that are required to call mutations. Please also check if `pon_sw.RDS` exists after running all steps shown under <workflow_PON>. Then run variant calling by mutect2 following the gwf workflow. 

***STEP4: Run VarScan2***

Running VarsCan2 requires mplieup file as input and doesnot require PON data. Please just follow the gwf workflow to perform VarScan2 variant calling.

***STEP5: Run DREAMS-vc***

Running DREAMS-vc requires the information extracted from PON and the model training by PON. Before running variant calling, please check whether `all_pon_soft.info.csv` and `all_pon_training_soft.vd.hdf5` exist after running all steps shown under <workflow_PON>. Then run variant calling by DREAMS-vc following the gwf workflow.


