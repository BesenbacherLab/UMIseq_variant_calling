import os
from gwf import AnonymousTarget
def fq_merge(fq_list, fq_dir, sample, options=None):
    '''
    Merge the fastq files (R1,R2,UMI) from different lanes for each sample
    '''
    inputs = [fq_list]
    outputs = [fq1,  fq2, fq3]

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
    mv ${{TEMP_DIR}}/tmp_R$((i+1)).fastq.gz {fq_dir}/{sample}/{sample}_R$((i+1)).fastq.gz
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


def tag_equalize(umicorrected_bam, umitagged_bam, options=None):
    '''
    From indexed bam, copy tag BX(corrected UMI), UG(unique id) from one read to the other read, add tag MQ(mapping quality)
    '''
    inputs = [umicorrected_bam]
    outputs = [umitagged_bam]
    if not options:
        options = dict(cores='16', memory='16g', walltime='4:00:00', account = 'ctDNA_var_calling')
    spec=f"""
    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}} 
    mkdir -p ${{TEMP_DIR}}

    samtools sort -n -T ${{TEMP_DIR}} -@ {options['cores']} -l 3 -o ${{TEMP_DIR}}/uminsorted.bam {umicorrected_bam}  
    python equalize_tags.py -q -t UG -t BX -i ${{TEMP_DIR}}/uminsorted.bam -o ${{TEMP_DIR}}/umitagged.bam

    mv ${{TEMP_DIR}}/umitagged.bam {umitagged_bam}
    """
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def fgbio_group(umitagged_bam, umigrouped_bam, min_umi_length, grouping_strategy, family_size_histogram, options=None):
    '''
    Group reads by the corrected umi made by umi_tools grouping
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
    '''
    Make a consensus read for each group resulted from umitools and filter reads 
    '''
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


def pe_ubam_mapper(input_ubam_file, output_bam_file, output_bai_file, reference_path, read_group, options=None):
    '''
    Remap the unmapped BAM file generated from fgbio consensus step, and copy the tags from unmapped BAM to the new BAM
    '''

    if not options:
        options = dict(cores='16', memory='16g',
                       walltime='12:00:00', account='ctdna_var_calling')

    inputs = [input_ubam_file]

    outputs = [output_bam_file, output_bai_file]

    threads = options['cores']
    rg = "@RG"
    for k, v in read_group.items():
        rg = rg + f"\\t{k}:{v}"
    spec = f'''

    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}}
    mkdir -p ${{TEMP_DIR}}

    gatk --java-options '-Xmx{options['memory']}' SortSam --TMP_DIR ${{TEMP_DIR}} -I {input_ubam_file} -SO queryname -OUTPUT ${{TEMP_DIR}}/unmapped_sorted.sam

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
    --OUTPUT /dev/stdout  \
    --TMP_DIR ${{TEMP_DIR}} | \
    
    gatk --java-options '-Xmx{options['memory']}' MergeBamAlignment \
    --ALIGNED_BAM /dev/stdin \
    --UNMAPPED_BAM ${{TEMP_DIR}}/unmapped_sorted.sam \
    --OUTPUT ${{TEMP_DIR}}/out.bam \
    -R {reference_path} \
    --INCLUDE_SECONDARY_ALIGNMENTS false \
    --CREATE_INDEX true \
    --TMP_DIR ${{TEMP_DIR}}
     
    mv ${{TEMP_DIR}}/out.bam {output_bam_file}
    mv ${{TEMP_DIR}}/out.bai {output_bai_file}
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
    '''
    Soft clip the overlapped regions to avoid double counnting variants 
    '''
    inputs = [realigned_consensus_bam]
    outputs = [clipped_consensus_bam, clipped_consensus_bai, clip_metrics]
    if not options:
        options = dict(cores='1', memory='4g', walltime='2:00:00', account = 'ctdna_var_calling')
    spec=f"""
    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}}
    mkdir -p ${{TEMP_DIR}}
    samtools sort -n -@ {options['cores']} -T ${{TEMP_DIR}} {realigned_consensus_bam} > ${{TEMP_DIR}}/nsorted.bam
    fgbio -Xmx{options['memory']} --tmp-dir=${{TEMP_DIR}} ClipBam \
        -i ${{TEMP_DIR}}/nsorted.bam \
        -o ${{TEMP_DIR}}/clipped_consensus.bam \
        -c SoftWithMask \
        -r {ref} \
        -m ${{TEMP_DIR}}/clip.metrics \
        --clip-overlapping-reads=true

    mv ${{TEMP_DIR}}/clipped_consensus.bam {clipped_consensus_bam}
    mv ${{TEMP_DIR}}/clipped_consensus.bai {clipped_consensus_bai}
    mv ${{TEMP_DIR}}/clip.metrics {clip_metrics}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, spec=spec, options=options)


def target_subset(clipped_consensus_bam, panel_bed, target_consensus_bam, options=None):
    '''
    Only output alignments overlapping the input BED FILE
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
    '''
    After overlapping clipping and target region subset, 
    some reads can lost their paired read but still keep the old tags, 
    so check their tags and assign the right tags for downstream filtering 
    '''
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
    '''
    Generate a PILEUP file from BAM file
    '''
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


def shearwater_new(and_all_vcf, and_all_tempvcf, or_all_vcf, or_all_tempvcf, pon, pileup_file, panel_bed, sample, options=None):
    '''
    Adapted from the DEEPSNV shearwater algorithm to make it work for large-sized data
    all_vcf has set up the threshold of 0.05 thus only outputs the significant variants
    temp_vcf outputs all types of variants along the whole region
    Work for both AND & OR algorithm 
    '''
    inputs = [pileup_file]
    outputs = [and_all_tempvcf, or_all_tempvcf]
    if not options:
        options = dict(cores='1', memory='32g', walltime='4:00:00', account = 'ctdna_var_calling')
    spec = f"""
    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}}
    mkdir -p ${{TEMP_DIR}}

    source ~/miniconda3/bin/activate shearwater
    Rscript shearwater.R \
        --vcfall=${{TEMP_DIR}}/and.all.vcf \
        --tempvcfall=${{TEMP_DIR}}/and.all.tempvcf \
        --pon={pon} \
        --region={panel_bed} \
        --model='AND' \
        --sample={sample} \
        {pileup_file}
    
    Rscript shearwater.R \
        --vcfall=${{TEMP_DIR}}/or.all.vcf \
        --tempvcfall=${{TEMP_DIR}}/or.all.tempvcf \
        --pon={pon} \
        --region={panel_bed} \
        --model='OR' \
        --sample={sample} \
        {pileup_file}

    mv ${{TEMP_DIR}}/and.all.vcf {and_all_vcf}
    mv ${{TEMP_DIR}}/and.all.tempvcf {and_all_tempvcf}

    mv ${{TEMP_DIR}}/or.all.vcf {or_all_vcf}
    mv ${{TEMP_DIR}}/or.all.tempvcf {or_all_tempvcf}
    
    """
    return AnonymousTarget(inputs = inputs, outputs = outputs, spec = spec, options = options)


def mutect2_filtergermline(ref, sorted_target_consensus_bam, mutation_positions_vcf, germline_sorted_target_bam, normal, pon_vcf, vcf, vcf_stats, options=None):
    '''
    Mutect2 with the input of germline BAM, outputs all types of variants along the whole region 
    '''
    inputs = [ref, sorted_target_consensus_bam, germline_sorted_target_bam]
    outputs = [vcf, vcf_stats]
    if not options:
        options = dict(cores='16', memory='32g', walltime='12:00:00', account = 'ctdna_var_calling')

    spec = f"""
    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}}
    mkdir -p ${{TEMP_DIR}}

    gatk --java-options '-Xmx{options['memory']} -XX:ParallelGCThreads={options['cores']}' Mutect2 \
        --tmp-dir ${{TEMP_DIR}} \
        -R {ref} \
        -I {sorted_target_consensus_bam} \
        -I {germline_sorted_target_bam} \
        -alleles {mutation_positions_vcf} \
        -normal "{normal}" \
        --panel-of-normals {pon_vcf} \
        --max-mnp-distance 0 \
        -O ${{TEMP_DIR}}/variant.vcf

    mv ${{TEMP_DIR}}/variant.vcf {vcf}
    mv ${{TEMP_DIR}}/variant.vcf.stats {vcf_stats}
    """
    return AnonymousTarget(inputs = inputs, outputs = outputs, spec = spec, options = options)


def mutect2_split(vcf, split_vcf, options=None):
    '''
    Change the mutect2 vcf format to the easy-analyzed common vcf format
    '''
    inputs = [vcf]
    outputs = [split_vcf]
    if not options:
        options = dict(cores='8', memory='8g', walltime='4:00:00', account = 'ctdna_var_calling')
    spec = f"""
    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}}
    mkdir -p ${{TEMP_DIR}}

    bcftools norm -Ov -m-any --atomize {vcf} > ${{TEMP_DIR}}/variant.split.vcf
    mv ${{TEMP_DIR}}/variant.split.vcf {split_vcf}
    """
    return AnonymousTarget(inputs = inputs, outputs = outputs, spec = spec, options = options)


def samtools_mpileup(bam,ref, min_BQ, mpileup, options=None):
    '''
    Use SAMTOOLS to make a PILEUP file required by Varscan2 
    '''
    inputs = [bam]
    outputs = [mpileup]
    if not options:
        options = dict(cores='4', memory='4g', walltime='4:00:00', account = 'ctdna_var_calling')
    spec = f"""
    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}}
    mkdir -p ${{TEMP_DIR}}
    
    samtools mpileup -f {ref} --max-depth 1000000 --min-BQ {min_BQ} {bam} > ${{TEMP_DIR}}/mpileup
    mv ${{TEMP_DIR}}/mpileup {mpileup}

    """
    return AnonymousTarget(inputs = inputs, outputs = outputs, spec = spec, options = options)


def varscan2(output_basename, mpileup, germline_mpileup, options=None):
    '''
    Run Varscan2 for variant calling 
    '''
    inputs = [mpileup, germline_mpileup]
    outputs = []
    if not options:
        options = dict(cores='4', memory='4g', walltime='4:00:00', account = 'ctdna_var_calling')
    spec = f"""
    
    varscan somatic {germline_mpileup} {mpileup} {output_basename} --validation 1 --somatic-p-value 1
    """
    return AnonymousTarget(inputs = inputs, outputs = outputs, spec = spec, options = options)


def dreams_variantcalls(mutations_df, bam, ref, model, all_pon_info, calls, options=None):
    '''
    Run DREAMS-vc for variant calling 
    '''
    inputs = [bam]
    outputs = [calls]
    if not options:
        options = dict(cores='1', memory='500g', walltime='24:00:00', account = 'ctdna_var_calling')
    spec = f"""
    set -e
    TEMP_DIR=temp/scratch/${{SLURM_JOBID}}
    mkdir -p ${{TEMP_DIR}}

    source ~/miniconda3/bin/activate dreams_M
    Rscript dreams3.R \
        {ref}\
        {mutations_df} \
        {bam} \
        {model} \
        {all_pon_info} \
        ${{TEMP_DIR}}/calls.df

    mv ${{TEMP_DIR}}/calls.df {calls}
    """
    return AnonymousTarget(inputs = inputs, outputs = outputs, spec = spec, options = options)

    






