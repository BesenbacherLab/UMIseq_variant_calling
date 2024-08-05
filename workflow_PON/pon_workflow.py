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
from templates import dreams_readdata
from gwf import Workflow, AnonymousTarget

fq_dir = param['meta']['pon_directory']
ref = param['references']['reference']
panel_bed = param['references']['panel_bed']

dirs = os.listdir(fq_dir)
pons=[]
for dir in dirs:
    if dir.startswith('Donor'):
        pons.append(dir)

gwf = Workflow()
for sample in pons:
    gwf.target_from_template(name = f"fqmerge_{sample}", template = fq_merge(fq_dir=fq_dir, sample=sample))
    
    gwf.target_from_template(name = f"fq_mapping_{sample}", 
        template = fq_mapping(bam=os.path.join(fq_dir,sample,f"{sample}.bam"),
            read_group={"ID":sample, "SM":sample, "PL":"illumina", "CN":"moma"},
            fq1 = os.path.join(fq_dir, sample, f"{sample}_R1.fastq.gz"),
            fq2 = os.path.join(fq_dir, sample, f"{sample}_R2.fastq.gz"),
            ref = ref))

    gwf.target_from_template(name = f"bam_sort_index_{sample}",
        template = bam_sort_index(bam = os.path.join(fq_dir, sample, f"{sample}.bam"), 
            sorted_bam = os.path.join(fq_dir, sample, f"{sample}_sorted.bam"),
            indexed_bam = os.path.join(fq_dir,sample,f"{sample}_sorted.bam.bai")))

    gwf.target_from_template(name = f"umitools_correct_{sample}", 
        template = umitools_correct(method="directional",
            sorted_bam=os.path.join(fq_dir, sample, f"{sample}_sorted.bam"),
            umicorrected_bam=os.path.join(fq_dir,sample,f"{sample}_umicorrected.bam"),
            umicorrected_tsv=os.path.join(fq_dir, sample,f"{sample}_umicorrected.tsv"),
            umicorrected_log=os.path.join(fq_dir,sample,f"{sample}_umicorrected.log")))

    gwf.target_from_template(name = f"tag_equalize_{sample}",
        template = tag_equalize(umicorrected_bam=os.path.join(fq_dir, sample, f"{sample}_umicorrected.bam"),
            uminsorted_bam=os.path.join(fq_dir, sample, f"{sample}_uminsorted.bam"),
            umitagged_bam=os.path.join(fq_dir, sample, f"{sample}_umitagged.bam")))

    gwf.target_from_template(name = f"fgbio_group_{sample}",
        template = fgbio_group(umitagged_bam=os.path.join(fq_dir, sample, f"{sample}_umitagged.bam"),
            umigrouped_bam=os.path.join(fq_dir, sample, f"{sample}_umigrouped.bam"), 
            min_umi_length=9, 
            grouping_strategy="identity", 
            family_size_histogram=os.path.join(fq_dir, sample, f"{sample}.hist.txt")))

    gwf.target_from_template(name= f"call_consensus_read_filter_{sample}",
        template = call_consensus_read_filter(umigrouped_bam=os.path.join(fq_dir, sample, f"{sample}_umigrouped.bam"), 
            rejected_bam=os.path.join(fq_dir,sample,f"{sample}_rejected.bam"), 
            unfiltered_consensus_bam=os.path.join(fq_dir,sample, f"{sample}_unfiltered.consensus.bam"), 
            filtered_consensus_bam=os.path.join(fq_dir, sample, f"{sample}_filtered.consensus.bam"), 
            filtered_consensus_bai=os.path.join(fq_dir, sample, f"{sample}_filtered.consensus.bai"),
            prefix=sample, 
            ref=ref, 
            min_reads=1, 
            max_read_error_rate=0.01, 
            max_base_error_rate=0.1, 
            min_base_quality=15, 
            max_no_call_fraction=0.2, 
            min_mean_base_quality=20))

    gwf.target_from_template(name=f"ubamtobam_{sample}",
        template=pe_ubam_mapper(input_ubam_file=os.path.join(fq_dir, sample, f"{sample}_filtered.consensus.bam"), 
            output_bam_file=os.path.join(fq_dir, sample, f"{sample}_utag.mapped.consensus.bam"), 
            read_group ={"ID":"A", "SM":sample, "PL":"illumina", "CN":"moma"},
            reference_path=ref))


    gwf.target_from_template(name = f"overlap_clip_{sample}",
        template=overlap_clip(ref=ref,
            realigned_consensus_bam=os.path.join(fq_dir, sample, f"{sample}_utag.mapped.consensus.bam"),
            clipped_consensus_bam=os.path.join(fq_dir, sample, f"{sample}_softclipped.consensus.bam"),
            clipped_consensus_bai=os.path.join(fq_dir, sample, f"{sample}_softclipped.consensus.bai"),
            clip_metrics=os.path.join(fq_dir, sample, "softclip.metrics"))) 

    gwf.target_from_template(name = f"target_subset_{sample}", 
        template=target_subset(clipped_consensus_bam=os.path.join(fq_dir, sample, f"{sample}_softclipped.consensus.bam"), 
            panel_bed=panel_bed, 
            target_consensus_bam=os.path.join(fq_dir, sample, f"{sample}_soft.target.consensus.bam")))  

    gwf.target_from_template(name = f"fix_mate_{sample}",
        template = fix_mate(bam = os.path.join(fq_dir, sample, f"{sample}_soft.target.consensus.bam"),
            filtered_fixmate_bam = os.path.join(fq_dir, sample, f"{sample}_soft.filtered.fixmate.bam"), 
            indexed_bam = os.path.join(fq_dir, sample, f"{sample}_soft.filtered.fixmate.bam.bai"))) 

    gwf.target_from_template(name = f"pileup_{sample}",
        template = pileup(input_bam = os.path.join(fq_dir, sample, f"{sample}_soft.filtered.fixmate.bam"),
            pileup_file = os.path.join(fq_dir, sample, f"{sample}_pileup"),
            ref = ref,
            min_baseq = 1))

    

    gwf.target_from_template(name = f"read_pon_{sample}",
        template = dreams_readdata(
            pon_bam = os.path.join(fq_dir, sample,f"{sample}_soft.filtered.fixmate.bam"), 
            ref = ref, 
            pon_data = os.path.join(fq_dir, sample,f"{sample}_soft.data.csv"), 
            pon_info = os.path.join(fq_dir, sample,f"{sample}_soft.info.csv")))