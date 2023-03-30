#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 shave_polish_4.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment pipeline for VEctor 
# Date:                 2022.10.05
# Run:                  snakemake -s workflow/rules/shave_polish_4.smk --profile slurm/
# Latest modification:  2023.01.25
# Done:                 Added cluster-config, get_threads(),

###############################################################################
# CONFIGURATION FILES #
configfile: "config/config.yaml"
cluster_config: "slurm/config.yaml"

import os, sys
from snakemake.utils import min_version
min_version("5.18.0")

###############################################################################
# WILDCARDS #
SAMPLE, = glob_wildcards("resources/reads/{sample}_R1.fastq.gz")

###############################################################################
# RESOURCES #
TMPDIR = config["resources"]["tmpdir"] # Temporary directory

###############################################################################
# PARAMETERS #
REFPATH = config["path"]                            # Path to genomes references
REFERENCE = config["reference"]                     # Genome reference sequence, in fasta format

###############################################################################
# FUNCTIONS AND COMMANDS #

# def get_threads(rule, default):
#     """
#     retrieve cpus-per-task value from cluster_config file available for SLURM
#     if fail, return default value defined on each rule

#     Example:
# 	rule bwa_mapping:
# 	    threads: get_threads('bwa_mapping', 8)
#     """
#     if rule in cluster_config and 'threads' in cluster_config[rule]:
#         return int(cluster_config[rule]['threads'])
#     elif rule in cluster_config and "cpus-per-task" in cluster_config[rule]:
#         return int(cluster_config[rule]["cpus-per-task"])
#     elif "__default__" in cluster_config and "cpus-per-task" in cluster_config["__default__"]:
#         return int(cluster_config["__default__"]["cpus-per-task"])
#     elif '__default__' in cluster_config and 'threads' in cluster_config['__default__']:
#         return int(cluster_config['__default__']['threads'])
#     return default

# def get_mem(rule, default):
#     """
#     retrieve mem_mb value form cluster_config file available for SLURM
#     if fail, return default value defined on each rule
#     Example:
#         rule bwa_mapping:
#             resources: get_mem('bwa_mapping', 8)
#     """
#     if rule in cluster_config and "mem-per-cpu" in cluster_config[rule]:
#         return int(cluster_config[rule]["mem-per-cpu"])
#     if "__default__" in cluster_config and "mem-per-cpu" in cluster_config["__default__"]:
#         return int(cluster_config["__default__"]["mem-per-cpu"])
#     return default

############################## O N S T A R T ###################################
onstart:
    shell("mkdir -p Cluster_logs/")

################################## A L L #######################################
rule all:
    input:
        multiqc = "multiqc_report.html",
        quickcheck = "bad_bams.fofn",
        check = expand("results/00_Quality_Control/validatesamfile/{sample}_bwa_md_realigned_fixed_ValidateSam.txt", sample=SAMPLE),
        flagstat = expand("results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed_bam.flagstat.txt", sample=SAMPLE),
        idxstats = expand("results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed.idxstats.txt", sample=SAMPLE),
        stats = expand("results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed_stats.txt", sample=SAMPLE),        

############################### R U L E S #####################################
rule samtools_quickcheck:
    input: 
        bam = "results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam",
    output:
        check = "bad_bams.fofn",
    resources: cpus=1, mem_mb=4000, time_min=120    
    log:
        "results/11_Reports/samtools/quickcheck/{sample}_bwa_realigned_fixed_bam.log",
    shell:
        config["MODULES"]["SAMTOOLS"]+"""
        samtools quickcheck -v {input.bam} >> {output.check} && echo 'all OK' || echo 'some bam files failed check, see bad_bams.fofn &> {log}
        """    
        
###############################################################################
rule samtools_stats:
    # Aim: Collects statistics from BAM files
    # Use: samtools stats -r ref.fa input.bam
    message:
        "SamTools stats"
    resources: cpus=1, mem_mb=4000, time_min=120
    input:
        bam = "results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam",
        ref = REFPATH+REFERENCE,
    output:
        stats = "results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed_stats.txt",
    log:
        "results/11_Reports/samtools/{sample}_bwa_realigned_fixed_stats.log",
    shell:
        config["MODULES"]["SAMTOOLS"]+"""
        samtools stats --threads {threads} -r {input.ref} {input.bam} 1> {output.stats} 2> {log}
        """

###############################################################################
rule validate_sam:
    # Aim: Basic check for bam file validity, as interpreted by the Broad Institute.
    # Use: picard.jar ValidateSamFile \
    #      -I input.bam \
    #      - MODE SUMMARY
    message:
        "Picard ValidateSamFile"
    resources: cpus=1, mem_mb=4000, time_min=120
    input:
        bam = "results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam",
    output:
        check = "results/00_Quality_Control/validatesamfile/{sample}_bwa_md_realigned_fixed_ValidateSam.txt",
    log:
        "results/11_Reports/validatesamfiles/{sample}_bwa_realigned_fixed_validate_bam.log",
    shell:
        config["MODULES"]["PICARDTOOLS"]+"""
        picard ValidateSamFile -I {input.bam} -O {output.check} -M SUMMARY > {log} 2>&1 || true
        """

###############################################################################
rule samtools_idxstats:
    # Aim: samtools idxstats – reports alignment summary statistics
    #       Retrieve and print stats in the index file corresponding to the input file. Before calling idxstats, the input BAM file should be indexed by samtools index.
    #       If run on a SAM or CRAM file or an unindexed BAM file, this command will still produce the same summary statistics, but does so by reading through the entire file.
    #       This is far slower than using the BAM indices.
    #       The output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped read-segments and # unmapped read-segments.
    #       It is written to stdout. Note this may count reads multiple times if they are mapped more than once or in multiple fragments.
    resources: cpus=1, mem_mb=4000, time_min=120 
    input:
        bam="results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam",
        idx="results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bai",
    output:
        "results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed.idxstats.txt",
    log:
        "results/11_Reports/samtools/idxstats/{sample}_bwa_realigned_fixed_idxstats.log",
    params:
        extra="",  # optional params string
    shell:
        config["MODULES"]["SAMTOOLS"]+"""
        samtools idxstats {input.bam} > {output} &> {log}
        """

###############################################################################
rule samtools_flagstat:
    input:
        bam = "results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam",
    output:
        flagstat = "results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed_bam.flagstat.txt",
    resources: cpus=1, mem_mb=4000, time_min=120    
    log:
        "results/11_Reports/samtools/flagstat/{sample}_bwa_realigned_fixed_bam.log",
    params:
        extra="",  # optional params string
    shell:
        config["MODULES"]["SAMTOOLS"]+"""
        samtools flagstat {input.bam} > {output.flagstat} &> {log}
        """

###############################################################################
rule multiqc:
    input:
        expand("results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed_bam.flagstat.txt", sample=SAMPLE),
        expand("results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed.idxstats.txt", sample=SAMPLE),
        expand("results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed_stats.txt", sample=SAMPLE),
        expand("results/02_Mapping/{sample}_bwa_sorted-mark-dup_metrics.txt", sample=SAMPLE),
        "results/00_Quality_Control/fastqc/",
        "results/00_Quality_Control/fastq-screen/",
        "results/00_Quality_Control/trimmed_fastqc/",
    output:
        "multiqc_report.html"
    resources: cpus=4, mem_mb=8000, time_min=120
    log:
        "results/11_Reports/multiqc/multiqc.log"
    shell:
        config["MODULES"]["MULTIQC"]+"""
        multiqc {input} -o results/00_Quality_Control/MULTIQC/ -n {output} > {log} 2>&1 
        """
