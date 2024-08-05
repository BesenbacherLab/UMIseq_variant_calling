# UMIseq_variant_calling
## Patient plasma sample processing
**About workflow**

All files under the folder <workflow> are used to process patient plasma sample files. It includes the UMI processing step and variant calling by VarScan2, Mutect2, shearwater and DREAMS-vc.

**About param.json**

Inside param.json, the working directory, conda environments, gwf workflow, addtional scripts for running, the reference files, the patterns for FASRQ files, and the parameters for UMI grouping and UMI consensus have been stated. Please check whether all the reference files exist before you run the workflow, and change the file paths and working directory for your own use. 

**Specfic steps** 

***STEP1:*** 

Run the gwf workflow with the templates.py and workflow.py under workflow.  
