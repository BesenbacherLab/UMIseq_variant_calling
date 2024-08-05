#!/bin/sh
gatk GenomicsDBImport -R /home/yixinlin/ctdna_var_calling/Workspaces/yixinlin/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
      --genomicsdb-workspace-path pon_db \
      -L /home/yixinlin/ctdna_var_calling/Workspaces/yixinlin/NEW_METHOD_hg38_08feb2016_capture_targets.bed \
      -V Donor284/Donor284.vcf \
      -V Donor285/Donor285.vcf \
      -V Donor286/Donor286.vcf \
      -V Donor287/Donor287.vcf \
      -V Donor288/Donor288.vcf \
      -V Donor289/Donor289.vcf \
      -V Donor290/Donor290.vcf \
      -V Donor291/Donor291.vcf \
      -V Donor292/Donor292.vcf \
      -V Donor293/Donor293.vcf \
      -V Donor294/Donor294.vcf \
      -V Donor295/Donor295.vcf \
      -V Donor296/Donor296.vcf \
      -V Donor297/Donor297.vcf \
      -V Donor298/Donor298.vcf \
      -V Donor299/Donor299.vcf \
      -V Donor300/Donor300.vcf \
      -V Donor301/Donor301.vcf \
      -V Donor302/Donor302.vcf \
      -V Donor303/Donor303.vcf \
      -V Donor304/Donor304.vcf \
      -V Donor305/Donor305.vcf \
      -V Donor306/Donor306.vcf \
      -V Donor307/Donor307.vcf \
      -V Donor308/Donor308.vcf \
      -V Donor309/Donor309.vcf \
      -V Donor310/Donor310.vcf \
      -V Donor311/Donor311.vcf \
      -V Donor312/Donor312.vcf \
      -V Donor313/Donor313.vcf \
      -V Donor314/Donor314.vcf \
      -V Donor315/Donor315.vcf \
      -V Donor316/Donor316.vcf \
      -V Donor317/Donor317.vcf \
      -V Donor318/Donor318.vcf \
      -V Donor319/Donor319.vcf \
      -V Donor320/Donor320.vcf \
      -V Donor321/Donor321.vcf \
      -V Donor322/Donor322.vcf \
      -V Donor323/Donor323.vcf \
      -V Donor324/Donor324.vcf \
      -V Donor325/Donor325.vcf \
      -V Donor326/Donor326.vcf \
      -V Donor327/Donor327.vcf \
      -V Donor328/Donor328.vcf \
      -V Donor329/Donor329.vcf 

gatk CreateSomaticPanelOfNormals -R /home/yixinlin/ctdna_var_calling/Workspaces/yixinlin/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
   -V gendb://pon_db \
   -O pon.vcf.gz

####after templates_pon and workflow_pon, step 2 for create pon.vcf for mutect2 filtering