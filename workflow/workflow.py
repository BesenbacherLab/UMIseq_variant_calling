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
from templates import pileup
from templates import shearwater_new
from templates import mutect2_filtergermline
from templates import mutect2_split
from templates import samtools_mpileup
from templates import varscan2
from templates import dreams_variantcalls
from templates import fix_mate
from gwf import Workflow, AnonymousTarget

with open('param.json') as f:
    param = json.load(f)

fq_dir = param['meta']['ctdna_directory'] #change to your own directory
germ_dir = param['meta']['germline_directory'] #change to your own directory
sample_file = param['meta']['sample_list'] #change to your own sample list
ref = param['references']['reference']
panel_bed = param['references']['panel_bed'] #change to your own bed file
pon_sw = param['references']['pon_sw']
pon_mutect = param['references']['pon_mutect']
pon_info_dreams = param["references"]["pon_info_dreams"]
positions_mutect = param['references']['positions_mutect']
positions_dreams = param['references']['positions_dreams']
model_dreams = param['references']['model_dreams']
grouping_method = param['parameters']['grouping_method']
min_consensus_reads = param['parameters']['min_consensus_reads'] 



gwf = Workflow()
for line in open(sample_file,'r'):
    sample_normal_tumor = line.strip().split(' ')
    sample = sample_normal_tumor[0]
    normal = sample_normal_tumor[1]
    tumor = sample_normal_tumor[2]
    ### Merge the fastq files(R1, R2, UMI) from different lanes for the same sample
    gwf.target_from_template(name = f"fqmerge_{sample}", 
        template = fq_merge(fq_dir = fq_dir, 
            sample = sample,
            fq_list = os.path.join(fq_dir, sample, f"{sample}_fq.lst"), #If there is no existing fq_list, you could run python prepare_fq_list.py to create fq_list for all samples
            fq1 = os.path.join(fq_dir, sample, f"{sample}_R1.fastq.gz"),
            fq2 = os.path.join(fq_dir, sample, f"{sample}_R2.fastq.gz"),
            fq3 = os.path.join(fq_dir, sample, f"{sample}_R3.fastq.gz")))
    
    ### Map the fastq files to hg38, get the BAM file
    gwf.target_from_template(name = f"fq_mapping_{sample}", 
        template = fq_mapping(bam=os.path.join(fq_dir,sample,f"{sample}.bam"),
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
            output_bai_file=os.path.join(fq_dir, sample, f"{sample}_utag.mapped.consensus.bai"),
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

    ### Make a PILEUP file from BAM for germline buffycoat sample
    gwf.target_from_template(name = f"germline_pileup_{sample}",
        template = pileup(input_bam = os.path.join(germ_dir, sample, f"{sample}_germline.filtered.fixmate.bam"),
            pileup_file = os.path.join(germ_dir, sample, f"{sample}_pileup"),
            ref = ref,
            min_baseq = 1)) 
         
    ### Variant calling by SHEARWATER
    gwf.target_from_template(name = f"shearwater_new_{sample}",
        template = shearwater_new(
            and_all_vcf = os.path.join(fq_dir, sample, f"{sample}_and_all.vcf"),
            and_all_tempvcf = os.path.join(fq_dir, sample, f"{sample}_and_all.tempvcf"),
            or_all_vcf = os.path.join(fq_dir, sample, f"{sample}_or_all.vcf"), 
            or_all_tempvcf = os.path.join(fq_dir, sample, f"{sample}_or_all.tempvcf"),
            pon = pon_sw, 
            pileup_file = os.path.join(fq_dir, sample, f"{sample}_pileup"), 
            panel_bed = panel_bed, 
            sample=sample))

    ### Variant calling by MUTECT2
    gwf.target_from_template(name = f"mutect2_filtergerm_{sample}",
        template = mutect2_filtergermline(
            sorted_target_consensus_bam = os.path.join(fq_dir, sample, f"{sample}_soft.filtered.fixmate.bam"),
            vcf = os.path.join(fq_dir, sample, f"{sample}_filgerm.vcf"),
            vcf_stats = os.path.join(fq_dir, sample, f"{sample}_filgerm.vcf.stats"),
            germline_sorted_target_bam = os.path.join(germ_dir, sample, f"{sample}_germline.filtered.fixmate.bam"),
            normal = normal,
            mutation_positions_vcf = positions_mutect,
            pon_vcf = pon_mutect,
            ref = ref))

    ### Change to an easy-analyzed format for MUTECT2 VCF
    gwf.target_from_template(name = f"mutect2_germsplit_{sample}",
        template = mutect2_split(
            vcf = os.path.join(fq_dir, sample, f"{sample}_filgerm.vcf"),
            split_vcf = os.path.join(fq_dir, sample, f"{sample}_germ.split.vcf")))

    ### Make a plasma PILEUP from BAM using samtools for running VARSCAN2
    gwf.target_from_template(name = f"ctdna_mpileup_{sample}",
        template = samtools_mpileup(
            bam = os.path.join(fq_dir, sample, f"{sample}_soft.filtered.fixmate.bam"),
            mpileup = os.path.join(fq_dir, sample, f"{sample}_mpileup"),
            min_BQ = 1,
            ref = ref))

    ### Make a germline buffycoat PILEUP from BAM using samtools for running VARSCAN2
    gwf.target_from_template(name = f"germline_mpileup_{sample}",
        template = samtools_mpileup(
            bam = os.path.join(germ_dir, sample, f"{sample}_germline.filtered.fixmate.bam"),
            mpileup = os.path.join(germ_dir, sample, f"{sample}_mpileup"),
            min_BQ = 1, 
            ref = ref))

    ### Variant calling by VARSCAN2
    gwf.target_from_template(name = f"varscan2_{sample}",
        template = varscan2(
            mpileup = os.path.join(fq_dir, sample, f"{sample}_mpileup"),
            germline_mpileup = os.path.join(germ_dir, sample, f"{sample}_mpileup"),
            output_basename = os.path.join(fq_dir, sample, f"{sample}_varscan")))

    ### Variant calling by DREAMS-vc
    gwf.target_from_template(name=f"dreams_variantcalls_{sample}",
        template = dreams_variantcalls(
            mutations_df = positions_dreams, 
            bam = os.path.join(fq_dir, sample, f"{sample}_soft.filtered.fixmate.bam"), 
            ref = ref, 
            model = model_dreams, 
            all_pon_info = pon_info_dreams, 
            calls = os.path.join(fq_dir, sample, f"{sample}_dreams_variants.df")))


