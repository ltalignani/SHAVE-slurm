#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 shave_markdup.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment pipeline for VEctor 
# Date:                 2022.10.05
# Run:                  snakemake -s workflow/rules/shave_markdup.smk --prfile slurm/
# Latest modification:  2023.02.03
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
ALIGNER = config["aligner"]                         # Aligners ('bwa' or 'bowtie2')

REFPATH = config["path"]                            # Path to genomes references
REFERENCE = config["reference"]                     # Genome reference sequence, in fasta format
INDEX = config["index"]                             # Genome reference index, in .fai format
DICTIONARY = config["dictionary"]                  # Genome reference dictionary, made w/ picard CreateSequenceDictionary, in .dict format

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
        bam = expand("results/02_Mapping/{sample}_bwa_sorted-mark-dup.bam", sample=SAMPLE),
        metrics=expand("results/02_Mapping/{sample}_bwa_sorted-mark-dup_metrics.txt", sample=SAMPLE),

############################### R U L E S #####################################
rule convert_and_sort:
    message:
        "Samtools view conversion of sample sam in bam format and sorting by coordinates"
    input:
        "results/02_Mapping/{sample}_bwa-mapped.sam",
    output:
        bam = "results/02_Mapping/{sample}_bwa_sorted.bam",
    resources: cpus=2, mem_mb=4000, time_min=120
    shell:
        config["MODULES"]["SAMTOOLS"]+"""
            (samtools view --threads {threads.cpus} -bS {input} | 
            samtools sort -o {output} -T {sample}.temp ) &> {log}
        """

###############################################################################
rule mark_duplicates_spark:
    message:
        "Picard MarkDuplicatesSpark remove PCR duplicates"
    input:
        "results/02_Mapping/{sample}_bwa_sorted.bam",
    output:
        bam = temp("results/02_Mapping/{sample}_bwa_sorted-mark-dup.bam"),
        metrics="results/02_Mapping/{sample}_bwa_sorted-mark-dup_metrics.txt",
    benchmark:
        "benchmarks/markduplicatesspark/{sample}_bwa.tsv"
    log:
        "results/11_Reports/markduplicatesspark/{sample}_bwa_sorted-mark-dup.log",
    params:
        extra="--remove-sequencing-duplicates",  # optional
        #java_opts=,  # optional
        #spark_runner="",  # optional, local by default
        #spark_v1.19.1="",  # optional
        #spark_extra="", # optional
    resources: cpus=16, mem_mb=16000, time_min=600
    shell:
        config["MODULES"]["GATK4"]+"\n"+config["MODULES"]["JAVA8"]+"""
            gatk MarkDuplicatesSpark -I {input} -O {output.bam} -M {output.metrics} {params.extra} > {log} 2>&1
        """




