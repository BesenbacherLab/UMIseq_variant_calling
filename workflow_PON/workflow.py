import os
import json
from templates import fq_merge
from templates import fq_mapping 
from templates import bam_sort_index
from templates import umitools_correct
from templates import tag_equalize
from templates import fgbio_group
from templates import call_consensus_read_filter
from templates import pe_ubam_mapper
from templates import overlap_clip
from templates import target_subset
from templates import fix_mate
from templates import pileup
from templates import pon_mutect2
from templates import vcf_index
from templates import create_pon_combined_vcf
from templates import dreams_readdata
from templates import dreams_trainmodel
from templates import bed_to_vcf
from templates import pon_shearwater
from gwf import Workflow, AnonymousTarget

with open('param.json') as f:
    param = json.load(f)

fq_dir = param['meta']['pon_directory'] #working directory
pon_db = param['meta']['pon_db'] 
ref = param['references']['reference']
panel_bed = param['references']['panel_bed'] #change to your own bed file
grouping_method = param['parameters']['grouping_method']
min_consensus_reads = param['parameters']['min_consensus_reads'] 

#each sub-directory(directory named by sample name) under working directory stores fastq files for each PON sample
dirs = os.listdir(fq_dir) 
pons=[] 
for dir in dirs:
    if dir.startswith('Donor'):
        pons.append(dir)

all_pon_pileup = []
all_pon_data = []
all_pon_info = []
all_pon_vcf = []
string = ''

gwf = Workflow()
for sample in pons:
    #####################################STEP1: UMI processing
    ### Merge the fastq files(R1, R2, UMI) from different lanes for the same PON sample
    gwf.target_from_template(name = f"fqmerge_{sample}", template = fq_merge(fq_dir=fq_dir, sample=sample))
    
    ### Map the fastq files to hg38, get the BAM file
    gwf.target_from_template(name = f"fq_mapping_{sample}", 
        template = fq_mapping(bam=os.path.join(fq_dir, sample,f"{sample}.bam"),
            read_group={"ID":sample, "SM":sample, "PL":"illumina", "CN":"moma"},
            fq1 = os.path.join(fq_dir, sample, f"{sample}_R1.fastq.gz"),
            fq2 = os.path.join(fq_dir, sample, f"{sample}_R2.fastq.gz"),
            ref = ref))
    
    ### Sort and index the BAM file
    gwf.target_from_template(name = f"bam_sort_index_{sample}",
        template = bam_sort_index(bam = os.path.join(fq_dir, sample, f"{sample}.bam"), 
            sorted_bam = os.path.join(fq_dir, sample, f"{sample}_sorted.bam"),
            indexed_bam = os.path.join(fq_dir,sample,f"{sample}_sorted.bam.bai")))

    ### Perform Umitools Grouping to output BAM with corrected UMI tag BX based on chosen grouping method
    gwf.target_from_template(name = f"umitools_correct_{sample}", 
        template = umitools_correct(method=grouping_method,
            sorted_bam=os.path.join(fq_dir, sample, f"{sample}_sorted.bam"),
            umicorrected_bam=os.path.join(fq_dir,sample,f"{sample}_umicorrected.bam"),
            umicorrected_tsv=os.path.join(fq_dir, sample,f"{sample}_umicorrected.tsv"),
            umicorrected_log=os.path.join(fq_dir,sample,f"{sample}_umicorrected.log")))
    
    ### Copy the tag BX from one read to the paired read
    gwf.target_from_template(name = f"tag_equalize_{sample}",
        template = tag_equalize(umicorrected_bam=os.path.join(fq_dir, sample, f"{sample}_umicorrected.bam"),
            uminsorted_bam=os.path.join(fq_dir, sample, f"{sample}_uminsorted.bam"),
            umitagged_bam=os.path.join(fq_dir, sample, f"{sample}_umitagged.bam")))
    
    ### Perform fgbio GroupReadsByUMI to group reads with identical BX together and assign tag MI for downstream steps
    gwf.target_from_template(name = f"fgbio_group_{sample}",
        template = fgbio_group(umitagged_bam=os.path.join(fq_dir, sample, f"{sample}_umitagged.bam"),
            umigrouped_bam=os.path.join(fq_dir, sample, f"{sample}_umigrouped.bam"), 
            min_umi_length=9, 
            grouping_strategy="identity", 
            family_size_histogram=os.path.join(fq_dir, sample, f"{sample}.hist.txt")))

    ### Make a consensus read from each read group and filter based on selected standards
    gwf.target_from_template(name= f"call_consensus_read_filter_{sample}",
        template = call_consensus_read_filter(umigrouped_bam=os.path.join(fq_dir, sample, f"{sample}_umigrouped.bam"), 
            rejected_bam=os.path.join(fq_dir,sample,f"{sample}_rejected.bam"), 
            unfiltered_consensus_bam=os.path.join(fq_dir,sample, f"{sample}_unfiltered.consensus.bam"), 
            filtered_consensus_bam=os.path.join(fq_dir, sample, f"{sample}_filtered.consensus.bam"), 
            filtered_consensus_bai=os.path.join(fq_dir, sample, f"{sample}_filtered.consensus.bai"),
            prefix=sample, 
            ref=ref, 
            min_reads=min_consensus_reads, 
            max_read_error_rate=0.01, 
            max_base_error_rate=0.1, 
            min_base_quality=15, 
            max_no_call_fraction=0.2, 
            min_mean_base_quality=20))
    
    ### Realign the derived unmapped BAM to hg38 and copy useful tags from the unmapped BAM to the mapped BAM
    gwf.target_from_template(name=f"ubamtobam_{sample}",
        template=pe_ubam_mapper(input_ubam_file=os.path.join(fq_dir, sample, f"{sample}_filtered.consensus.bam"), 
            output_bam_file=os.path.join(fq_dir, sample, f"{sample}_utag.mapped.consensus.bam"), 
            read_group ={"ID":"A", "SM":sample, "PL":"illumina", "CN":"moma"},
            reference_path=ref))

    ### Softclip the overlap regions to avoid double counting of mutations
    gwf.target_from_template(name = f"overlap_clip_{sample}",
        template=overlap_clip(ref=ref,
            realigned_consensus_bam=os.path.join(fq_dir, sample, f"{sample}_utag.mapped.consensus.bam"),
            clipped_consensus_bam=os.path.join(fq_dir, sample, f"{sample}_softclipped.consensus.bam"),
            clipped_consensus_bai=os.path.join(fq_dir, sample, f"{sample}_softclipped.consensus.bai"),
            clip_metrics=os.path.join(fq_dir, sample, "softclip.metrics"))) 
    
    ### Subset the target region by your BED file
    gwf.target_from_template(name = f"target_subset_{sample}", 
        template=target_subset(clipped_consensus_bam=os.path.join(fq_dir, sample, f"{sample}_softclipped.consensus.bam"), 
            panel_bed=panel_bed, 
            target_consensus_bam=os.path.join(fq_dir, sample, f"{sample}_soft.target.consensus.bam")))  

    ### Recheck and correct the tags for filtering
    gwf.target_from_template(name = f"fix_mate_{sample}",
        template = fix_mate(bam = os.path.join(fq_dir, sample, f"{sample}_soft.target.consensus.bam"),
            filtered_fixmate_bam = os.path.join(fq_dir, sample, f"{sample}_soft.filtered.fixmate.bam"), 
            indexed_bam = os.path.join(fq_dir, sample, f"{sample}_soft.filtered.fixmate.bam.bai"))) 

    ### Make a PILEUP file from BAM for plasma sample
    gwf.target_from_template(name = f"pileup_{sample}",
        template = pileup(input_bam = os.path.join(fq_dir, sample, f"{sample}_soft.filtered.fixmate.bam"),
            pileup_file = os.path.join(fq_dir, sample, f"{sample}_pileup"),
            ref = ref,
            min_baseq = 1))
    all_pon_pileup.append(os.path.join(fq_dir, sample, f"{sample}_pileup"))



########################################STEP2: Prepare for variant calling
    ### For DREAMS-vc, get pon data info
    gwf.target_from_template(name = f"read_pon_{sample}",
        template = dreams_readdata(
            pon_bam = os.path.join(fq_dir, sample,f"{sample}_soft.filtered.fixmate.bam"), 
            ref = ref, 
            pon_data = os.path.join(fq_dir, sample,f"{sample}_soft.data.csv"), 
            pon_info = os.path.join(fq_dir, sample,f"{sample}_soft.info.csv")))
    all_pon_data.append(os.path.join(fq_dir, sample,f"{sample}_soft.data.csv"))
    all_pon_data.append(os.path.join(fq_dir, sample,f"{sample}_soft.info.csv"))
    
    ### For Mutect2, create pon vcf and index
    gwf.target_from_template(name = f"pon_mutect_vcf_{sample}",
        template = pon_mutect2(
            pon_bam = os.path.join(fq_dir, sample, f"{sample}_soft.filtered.fixmate.bam"),
            pon_vcf = os.path.join(fq_dir, sample, f"{sample}.vcf"),  
            ref = ref))
    
    gwf.target_from_template(name = f"pon_vcf_index_{sample}",
        template = vcf_index(
            pon_vcf = os.path.join(fq_dir, sample, f"{sample}.vcf"), 
            pon_vcf_idx = os.path.join(fq_dir, sample, f"{sample}.vcf.idx")))
    string = string + "-V " + f"{fq_dir}/{sample}/{sample}.vcf "
    all_pon_vcf.append(os.path.join(fq_dir, sample, f"{sample}.vcf"))


###For mutect2, create combined pon.vcf.gz for mutect2 variant calling
gwf.target_from_template(name = "create_pon_combined_vcf",
    template = create_pon_combined_vcf(
        panel_bed = panel_bed, 
        ref = ref,
        pon_db = pon_db,
        all_pon_vcf = all_pon_vcf,
        string = string,
        pon_combined_vcf = os.path.join(fq_dir, "pon.vcf.gz")))

### For DREAMS-vc, train pon data
gwf.target_from_template(name = "train_pon",
    template = dreams_trainmodel(
        all_pon_data = all_pon_data, 
        all_pon_info = all_pon_info,
        model_file = os.path.join(fq_dir, 'all_pon_training_soft.vd.hdf5'),
        log_file = os.path.join(fq_dir, 'all_pon_training_soft.vd.log')))

### For mutect2 and dreams, convert bed into vcf for variant calling on all alts on all positions
gwf.target_from_template(name = "bed_to_vcf",
    template = bed_to_vcf( 
        panel_bed = panel_bed,
        pileup = os.path.join(fq_dir, pons[0], f"{pons[0]}_pileup"), #any existing pileup to input REF
        mutect_allpos_vcf = os.path.join(fq_dir, 'mutect_allpos.vcf'),
        mutect_allpos_vcf_idx = os.path.join(fq_dir, 'mutect_allpos.vcf.idx'),
        dreams_allpos_vcf = os.path.join(fq_dir, 'dreams_allpon.vcf')))

### For shearwater, save all pon pileup counts into pon_sw.RDS for shearwater variant calling
gwf.target_from_template(name = "pon_sw",
    template = pon_shearwater(
        all_pon_pileup = all_pon_pileup, 
        panel_bed = panel_bed,
        pon_counts = os.path.join(fq_dir, 'pon_sw.RDS')))
