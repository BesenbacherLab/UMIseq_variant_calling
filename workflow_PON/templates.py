import os
from gwf import AnonymousTarget
def fq_merge(fq_dir, sample, options=None):
    '''
    Merge the fastq files (R1,R2,UMI) from different lanes for each sample
    '''
    fq_list = os.path.join(fq_dir, sample, f"{sample}_fq.lst")
    inputs = [fq_list]
    outputs = [os.path.join(fq_dir, sample, f"{sample}_R1.fastq.gz"),  os.path.join(fq_dir, sample, f"{sample}_R2.fastq.gz"), os.path.join(fq_dir, sample, f"{sample}_R3.fastq.gz")]

    if not options:
        options = dict(cores = '1', memory = '8g', walltime = '2:00:00', account = 'ctdna_var_calling')

    spec=f"""
    set -e 
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}} 
    mkdir -p ${{TEMP_DIR}}

    while read -a line; do 
        for i in ${{!line[@]}}; do 
            cat ${{line[i]}} >> ${{TEMP_DIR}}/tmp_R$((i+1)).fastq.gz 
        done 
    done < {fq_list}
    mv ${{TEMP_DIR}}/tmp_R1.fastq.gz {fq_dir}/{sample}/{sample}_R1.fastq.gz
    mv ${{TEMP_DIR}}/tmp_R2.fastq.gz {fq_dir}/{sample}/{sample}_R2.fastq.gz
    mv ${{TEMP_DIR}}/tmp_R3.fastq.gz {fq_dir}/{sample}/{sample}_R3.fastq.gz
    touch {fq_dir}/{sample}/fqmerge.done
    """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def fq_mapping(bam, read_group, fq1, fq2, ref, options=None):
    '''
    Map two paired fastq files against hg38 and get corresponding bam file
    '''
    inputs = [fq1, fq2]

    outputs = [bam, bam+".md5"] 
    if not options:
        options = dict(cores='16', memory='16g', walltime='6:00:00',  account = 'ctdna_var_calling')
    rg = "@RG"
    for k, v in read_group.items():
        rg = rg + f"\\t{k}:{v}"
    spec = f"""
    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}} 
    mkdir -p ${{TEMP_DIR}}

    bwa mem -K 100000000 -R '{rg}' -M -t {options["cores"]} {ref} {fq1} {fq2} \
    | samtools view -F 256 -F 4 -F 8 -q 1 -bo ${{TEMP_DIR}}/temp.bam
    
    samtools flagstat ${{TEMP_DIR}}/temp.bam | md5sum | cut -f1 -d " " > ${{TEMP_DIR}}/temp.bam.md5

    mv ${{TEMP_DIR}}/temp.bam {bam}
    mv ${{TEMP_DIR}}/temp.bam.md5 {bam}.md5
    """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def umitools_correct(method, sorted_bam, umicorrected_bam, umicorrected_tsv, umicorrected_log, options=None):
    '''
    Add real UMI tag(BX) for each read in a sample, output a tagged bam and a specific explanation(tsv) 
    '''
    inputs = [sorted_bam]
    outputs = [umicorrected_bam, umicorrected_tsv, umicorrected_log, umicorrected_bam + '.md5']
       
    if not options:
        options = dict(cores='1', memory='30g', walltime='4:00:00', account = 'ctDNA_var_calling')
    spec = f"""
    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}}
    mkdir -p ${{TEMP_DIR}}
    
    umi_tools group \
        --temp-dir=${{TEMP_DIR}} \
        --output-bam \
        --extract-umi-method=read_id \
        --umi-separator=":" \
        --umi-group-tag=BX \
        --method={method} \
        --edit-distance-threshold=1 \
        --paired \
        --chimeric-pairs=discard \
        --unpaired-reads=discard \
        -L ${{TEMP_DIR}}/umicorrected.log \
        -I {sorted_bam} \
        --group-out=${{TEMP_DIR}}/umicorrected.tsv \
        -S ${{TEMP_DIR}}/umicorrected.bam  

    samtools flagstat ${{TEMP_DIR}}/umicorrected.bam | md5sum | cut -f1 -d " " > ${{TEMP_DIR}}/umicorrected.bam.md5

    mv ${{TEMP_DIR}}/umicorrected.log {umicorrected_log}
    mv ${{TEMP_DIR}}/umicorrected.bam {umicorrected_bam}
    mv ${{TEMP_DIR}}/umicorrected.tsv {umicorrected_tsv}
    mv ${{TEMP_DIR}}/umicorrected.bam.md5 {umicorrected_bam}.md5
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, spec=spec, options=options)



def tag_equalize(umicorrected_bam, uminsorted_bam, umitagged_bam, options=None):
    '''
    From indexed bam, copy tag BX(corrected UMI), UG(unique id) from one read to the other read, add tag MQ(mapping quality)
    '''
    inputs = [umicorrected_bam]
    outputs = [umitagged_bam, uminsorted_bam]
    if not options:
        options = dict(cores='16', memory='16g', walltime='4:00:00', account = 'ctDNA_var_calling')
    spec=f"""
    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}} 
    mkdir -p ${{TEMP_DIR}}

    samtools sort -n -T ${{TEMP_DIR}} -l 3 -o ${{TEMP_DIR}}/uminsorted.bam {umicorrected_bam}  
    python equalize_tags.py -q -t UG -t BX -i ${{TEMP_DIR}}/uminsorted.bam -o ${{TEMP_DIR}}/umitagged.bam

    mv ${{TEMP_DIR}}/uminsorted.bam {uminsorted_bam}
    mv ${{TEMP_DIR}}/umitagged.bam {umitagged_bam}
    """
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def fgbio_group(umitagged_bam, umigrouped_bam, min_umi_length, grouping_strategy, family_size_histogram, options=None):
    '''
    group reads by the corrected umi from umi_tools group performance
    set -e: Exit immediately if a command exits with a non-zero
    mkdir -p: create the directory and, if required, all parent directories
    '''
    inputs = [umitagged_bam]
    outputs = [umigrouped_bam, family_size_histogram]
    if not options:
        options = dict(cores='4', memory='32g', walltime='8:00:00', account = 'ctdna_var_calling')
    spec=f"""
    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}} 
    mkdir -p ${{TEMP_DIR}}
    fgbio -Xmx{options['memory']} --tmp-dir=${{TEMP_DIR}} GroupReadsByUmi \
        -i {umitagged_bam} \
        -o ${{TEMP_DIR}}/umigrouped.bam \
        -s {grouping_strategy} \
        -f ${{TEMP_DIR}}/family_size_histogram \
        -t BX \
        -T MI \
        -n False \
        -e 0 \
        -l {min_umi_length}
    mv ${{TEMP_DIR}}/family_size_histogram {family_size_histogram}
    mv ${{TEMP_DIR}}/umigrouped.bam {umigrouped_bam}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def call_consensus_read_filter(umigrouped_bam, rejected_bam, unfiltered_consensus_bam, filtered_consensus_bam, filtered_consensus_bai, prefix, ref, min_reads, max_read_error_rate, max_base_error_rate, min_base_quality, max_no_call_fraction, min_mean_base_quality, options=None):
    inputs = [umigrouped_bam]
    outputs = [rejected_bam, unfiltered_consensus_bam, filtered_consensus_bai, filtered_consensus_bam, filtered_consensus_bam + '.md5']
    if not options:
        options = dict(cores='4', memory='32g', walltime='8:00:00', account = 'ctdna_var_calling')

    spec = f"""
    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}} 
    mkdir -p ${{TEMP_DIR}}
    fgbio -Xmx{options['memory']} --tmp-dir=${{TEMP_DIR}} CallMolecularConsensusReads \
        -t MI \
        -i {umigrouped_bam}\
        -o ${{TEMP_DIR}}/unfiltered_consensus.bam \
        -r ${{TEMP_DIR}}/rejected.bam \
        -p {prefix} \
        --error-rate-pre-umi=45 \
        --error-rate-post-umi=40 \
        --min-input-base-quality=10 \
        --min-reads=1 \
        --max-reads=1000 \
        --output-per-base-tags=true
    
    fgbio -Xmx{options['memory']} --tmp-dir=${{TEMP_DIR}} FilterConsensusReads \
          -i ${{TEMP_DIR}}/unfiltered_consensus.bam \
          -o ${{TEMP_DIR}}/filtered_consensus.bam \
          -r {ref} \
          --reverse-per-base-tags=true \
          --min-reads={min_reads} \
          --max-base-error-rate={max_base_error_rate} \
          --max-read-error-rate={max_read_error_rate} \
          --min-base-quality={min_base_quality} \
          --min-mean-base-quality={min_mean_base_quality} \
          --max-no-call-fraction={max_no_call_fraction}

    samtools flagstat ${{TEMP_DIR}}/filtered_consensus.bam | md5sum | cut -f1 -d " " > ${{TEMP_DIR}}/filtered_consensus.bam.md5

    mv ${{TEMP_DIR}}/unfiltered_consensus.bam {unfiltered_consensus_bam}
    mv ${{TEMP_DIR}}/rejected.bam {rejected_bam}
    mv ${{TEMP_DIR}}/filtered_consensus.bam {filtered_consensus_bam}
    mv ${{TEMP_DIR}}/filtered_consensus.bai {filtered_consensus_bai}
    mv ${{TEMP_DIR}}/filtered_consensus.bam.md5 {filtered_consensus_bam}.md5
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, spec=spec, options=options)


def pe_ubam_mapper(input_ubam_file, output_bam_file, reference_path, read_group, options=None):
    '''
    Map paired ends reads from umapped bam file. Suffix is used as read group ID.
    '''

    if not options:
        options = dict(cores='16', memory='16g',
                       walltime='6:00:00', account='ctdna_var_calling')

    inputs = [input_ubam_file]

    outputs = [output_bam_file]

    threads = options['cores']
    rg = "@RG"
    for k, v in read_group.items():
        rg = rg + f"\\t{k}:{v}"
    spec = f'''

    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}}
    mkdir -p ${{TEMP_DIR}}

    gatk --java-options '-Xmx{options['memory']}' SortSam -I {input_ubam_file} -SO queryname -OUTPUT ${{TEMP_DIR}}/unmapped_sorted.sam

    gatk --java-options '-Xmx{options['memory']}' SamToFastq \
    -I ${{TEMP_DIR}}/unmapped_sorted.sam \
    --FASTQ /dev/stdout \
    --INTERLEAVE true \
    --NON_PF true \
    --TMP_DIR ${{TEMP_DIR}} | \
    
    bwa mem -M -R "{rg}" -t {threads} -p {reference_path} /dev/stdin | \
    
    samtools view -hb -F 256 | \
    
    gatk --java-options '-Xmx{options['memory']}' SortSam \
    -I /dev/stdin \
    -SO queryname \
    --OUTPUT /dev/stdout | \
    
    gatk --java-options '-Xmx{options['memory']}' MergeBamAlignment \
    --ALIGNED_BAM /dev/stdin \
    --UNMAPPED_BAM ${{TEMP_DIR}}/unmapped_sorted.sam \
    --OUTPUT ${{TEMP_DIR}}/out.bam \
    -R {reference_path} \
    --INCLUDE_SECONDARY_ALIGNMENTS false \
    --CREATE_INDEX true \
    --TMP_DIR ${{TEMP_DIR}}
     
    mv ${{TEMP_DIR}}/out.bam {output_bam_file}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, spec=spec, options=options)


def bam_sort_index(bam, sorted_bam, indexed_bam, options=None):
    '''
    Prepare sorted and indexed bam file 
    '''
    inputs = [bam]
    outputs = [sorted_bam, indexed_bam]
    if not options:
        options = dict(cores='16', memory='16g', walltime='4:00:00', account = 'ctDNA_var_calling')
    spec = f"""
    set -e 
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}}
    mkdir -p ${{TEMP_DIR}}
    samtools sort -@ {options['cores']} -T ${{TEMP_DIR}} -o ${{TEMP_DIR}}/sorted.bam {bam}
    samtools index -@ {options['cores']} ${{TEMP_DIR}}/sorted.bam

    mv ${{TEMP_DIR}}/sorted.bam {sorted_bam}
    mv ${{TEMP_DIR}}/sorted.bam.bai {indexed_bam}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, spec=spec, options=options)


def overlap_clip(ref, realigned_consensus_bam, clipped_consensus_bam, clipped_consensus_bai, clip_metrics, options=None):
    inputs = [realigned_consensus_bam]
    outputs = [clipped_consensus_bam, clipped_consensus_bai, clip_metrics]
    if not options:
        options = dict(cores='1', memory='4g', walltime='2:00:00', account = 'ctdna_var_calling')
    spec=f"""
    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}}
    mkdir -p ${{TEMP_DIR}}
    fgbio -Xmx{options['memory']} --tmp-dir=${{TEMP_DIR}} ClipBam \
        -i {realigned_consensus_bam} \
        -o ${{TEMP_DIR}}/clipped_consensus.bam \
        -r {ref} \
        -c SoftWithMask \
        -m ${{TEMP_DIR}}/clip.metrics \
        --clip-overlapping-reads=true

    mv ${{TEMP_DIR}}/clipped_consensus.bam {clipped_consensus_bam}
    mv ${{TEMP_DIR}}/clipped_consensus.bai {clipped_consensus_bai}
    mv ${{TEMP_DIR}}/clip.metrics {clip_metrics}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, spec=spec, options=options)


def target_subset(clipped_consensus_bam, panel_bed, target_consensus_bam, options=None):
    '''
    -L: Only output alignments overlapping the input BED FILE
    '''
    inputs = [clipped_consensus_bam]
    outputs = [target_consensus_bam]
    if not options:
        options = dict(cores='8', memory='16g', walltime='12:00:00', account = 'ctdna_var_calling')
    
    spec=f"""
    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}}
    mkdir -p ${{TEMP_DIR}}
    samtools view {clipped_consensus_bam} -L {panel_bed} -o ${{TEMP_DIR}}/target_consensus.bam

    mv ${{TEMP_DIR}}/target_consensus.bam {target_consensus_bam}
    """
    return AnonymousTarget(inputs = inputs, outputs=outputs, options=options, spec=spec)


def fix_mate(bam, filtered_fixmate_bam, indexed_bam, options=None):
    inputs = [bam]
    outputs = [filtered_fixmate_bam, indexed_bam]
    if not options:
        options = dict(cores='4', memory='12g', walltime='4:00:00', account = 'ctdna_var_calling')
    spec = f"""
    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}}
    mkdir -p ${{TEMP_DIR}}

    samtools sort -T ${{TEMP_DIR}} -@ {options['cores']} -n -o ${{TEMP_DIR}}/nsort.bam {bam}
    samtools fixmate -@ {options['cores']} ${{TEMP_DIR}}/nsort.bam ${{TEMP_DIR}}/fixmate.bam
    samtools sort -T ${{TEMP_DIR}} -@ {options['cores']} ${{TEMP_DIR}}/fixmate.bam | samtools view -f 1 -F 780 -bho ${{TEMP_DIR}}/filtered.fixmate.bam 
    
    samtools index -@ {options['cores']} ${{TEMP_DIR}}/filtered.fixmate.bam

    mv ${{TEMP_DIR}}/filtered.fixmate.bam {filtered_fixmate_bam}
    mv ${{TEMP_DIR}}/filtered.fixmate.bam.bai {indexed_bam}
    """
    return AnonymousTarget(inputs = inputs, outputs = outputs, spec = spec, options = options)



def pileup(ref, input_bam, pileup_file, min_baseq, options=None):
    inputs = [input_bam]
    outputs = [pileup_file, pileup_file+".md5"]
    if not options:
        options = dict(cores='1', memory='4g', walltime='4:00:00', account ='ctdna_var_calling')
    spec = f"""
    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}}
    mkdir -p ${{TEMP_DIR}}
    
    pysamstats -f {ref} --output=${{TEMP_DIR}}/pileup --type=variation_strand \
        --max-depth=1000000 --min-baseq={min_baseq} {input_bam}

    md5sum ${{TEMP_DIR}}/pileup > ${{TEMP_DIR}}/pileup.md5

    mv ${{TEMP_DIR}}/pileup {pileup_file}
    mv ${{TEMP_DIR}}/pileup.md5 {pileup_file}.md5
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def pon_mutect2(ref, pon_bam, pon_vcf, options=None):
    inputs = [pon_bam]
    outputs = [pon_vcf]
    if not options:
        options = dict(cores='16', memory='32g', walltime='4:00:00', account = 'ctdna_var_calling')
    spec = f"""
    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}}
    mkdir -p ${{TEMP_DIR}}

    gatk --java-options '-Xmx{options['memory']}' Mutect2 \
        --tmp-dir ${{TEMP_DIR}} \
        -R {ref} \
        -I {pon_bam} \
        -max-mnp-distance 0 \
        -O ${{TEMP_DIR}}/pon.vcf

    mv ${{TEMP_DIR}}/pon.vcf {pon_vcf}
    """
    return AnonymousTarget(inputs = inputs, outputs = outputs, spec = spec, options = options)

def vcf_index(pon_vcf, pon_vcf_idx, options=None):
    inputs = [pon_vcf]
    outputs = [pon_vcf_idx]
    if not options:
        options = dict(cores='16', memory='32g', walltime='4:00:00', account = 'ctdna_var_calling')
    spec = f"""
    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}}
    mkdir -p ${{TEMP_DIR}}

    gatk --java-options '-Xmx{options['memory']}' IndexFeatureFile \
        --tmp-dir ${{TEMP_DIR}} \
        -I {pon_vcf} \
        -O ${{TEMP_DIR}}/pon.vcf.idx

    mv ${{TEMP_DIR}}/pon.vcf.idx {pon_vcf_idx}
    """
    return AnonymousTarget(inputs = inputs, outputs = outputs, spec = spec, options = options)


def dreams_readdata(pon_bam, ref, pon_data, pon_info, options=None):
    inputs = [pon_bam]
    outputs = [pon_data, pon_info]
    if not options:
        options = dict(cores='1', memory='25g', walltime='4:00:00', account = 'ctdna_var_calling')
    spec = f"""
    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}}
    mkdir -p ${{TEMP_DIR}}

    source ~/miniconda3/bin/activate dreams_M
    Rscript dreams.R \
        --pondata=${{TEMP_DIR}}/pon_data.csv \
        --poninfo=${{TEMP_DIR}}/pon_info.csv \
        --ref={ref} \
        {pon_bam}

    mv ${{TEMP_DIR}}/pon_data.csv {pon_data}
    mv ${{TEMP_DIR}}/pon_info.csv {pon_info}
    """
    return AnonymousTarget(inputs = inputs, outputs = outputs, spec = spec, options = options)