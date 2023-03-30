#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 shave_polish_1.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment pipeline for VEctor 
# Date:                 2022.10.05
# Run:                  snakemake -s workflow/rules/shave_polish_1.smk --profile slurm/
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
#     retrieve mem-per-cpu value form cluster_config file available for SLURM
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
        fix = expand("results/02_Mapping/{sample}_bwa_sorted-mark-dup-fx.bam", sample=SAMPLE),
        index = expand("results/02_Mapping/{sample}_bwa_sorted-mark-dup-fx.bai", sample=SAMPLE),

############################### R U L E S #####################################
rule samtools_index_markdup:
    # Aim: indexing marked as duplicate BAM file
    # Use: samtools index -@ [THREADS] -b [MARKDUP.bam] [INDEX.bai]
    message:
        "SamTools indexing marked as duplicate BAM file"
    resources: cpus=2, mem_mb=4000, time_min=120
    input:
        markdup = "results/02_Mapping/{sample}_bwa_sorted-mark-dup-fx.bam",
    output:
        index = "results/02_Mapping/{sample}_bwa_sorted-mark-dup-fx.bai",
    log:
        "results/11_Reports/samtools/{sample}_bwa_sorted-mark-dup-index.log",
    shell:
        config["MODULES"]["SAMTOOLS"]+"""
            samtools index -@ {threads} -b {input.markdup} {output.index} &> {log}
        """

###############################################################################
rule SetNmMdAndUqTags:
    # Aim: This tool takes in a coordinate-sorted SAM or BAM and calculates the NM, MD, and UQ tags by comparing with the reference.
    # Use: picard.jar SetNmMdAndUqTags \
    #       R=reference_sequence.fasta
    #       I=sorted.bam \
    #       O=fixed.bam
    message:
        "Picard SetNmMdAndUqTags"
    input:
        bam = "results/02_Mapping/{sample}_bwa_sorted-mark-dup.bam",
        ref = REFPATH+REFERENCE,
    output:
        fix = "results/02_Mapping/{sample}_bwa_sorted-mark-dup-fx.bam",
    resources: cpus=8, mem_mb=4000, time_min=120
    log:
        "results/11_Reports/SetNmMdAndUqTags/{sample}_bwa_sorted-mark-dup-fx.log",
    benchmark:
        "benchmarks/setnmmdanduqtags/{sample}_bwa.tsv",
    shell:
        config["MODULES"]["PICARDTOOLS"]+"""
            picard SetNmMdAndUqTags R={input.ref} I={input.bam} O={output.fix} > {log} 2>&1
        """