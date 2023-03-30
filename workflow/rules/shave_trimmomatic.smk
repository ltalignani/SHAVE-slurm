#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 shave_trimmomatic.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment pipeline for VEctor 
# Date:                 2022.10.05
# Run:                  snakemake -s workflow/rules/shave_trimmomatic.smk --profile slurm/
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
CONFIG = config["fastq-screen"]["config"]           # Fastq-screen --conf
MAPPER = config["fastq-screen"]["aligner"]          # Fastq-screen --aligner
SUBSET = config["fastq-screen"]["subset"]           # Fastq-screen --subset

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
#     retrieve mem-per-cpus value form cluster_config file available for SLURM
#     if fail, return default value defined on each rule
#     Example:
#         rule bwa_mapping:
#             resources: get_mem('bwa_mapping', 8)
#     """
#     if rule in cluster_config and "mem-per-cpus" in cluster_config[rule]:
#         return int(cluster_config[rule]["mem-per-cpus"])
#     if "__default__" in cluster_config and "mem-per-cpus" in cluster_config["__default__"]:
#         return int(cluster_config["__default__"]["mem-per-cpus"])
#     return default

############################## O N S T A R T ###################################

################################## A L L #######################################
rule all:
    input:
        fastqc          =   "results/00_Quality_Control/fastqc/", #expand("results/00_Quality_Control/fastqc/{sample}_fastqc.html", sample=SAMPLE),
        trimmed_fastqc  =   "results/00_Quality_Control/trimmed_fastqc/",
        fastqscreen     =   "results/00_Quality_Control/fastq-screen/",
        forward_reads   =   expand("results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R1.fastq.gz", sample=SAMPLE),
        reverse_reads   =   expand("results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R2.fastq.gz", sample=SAMPLE),
        forwardUnpaired =   expand("results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_unpaired_R1.fastq.gz", sample=SAMPLE),
        reverseUnpaired =   expand("results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_unpaired_R2.fastq.gz", sample=SAMPLE),

############################### R U L E S #####################################
rule fastqc_quality_control:
    #threads: get_threads('fastqc_quality_control', 1)
    message: "FastQC reads quality controling"
    resources: cpus=1, mem_mb=4000, tim_min=60                              #get_mem('fastqc_quality_control', 4000)
    input:
        fastq = "resources/reads/"
    output:
        fastqc = directory("results/00_Quality_Control/fastqc/")
    log:
        "results/11_Reports/quality/fastqc.log"
    shell:
        config["MODULES"]["FASTQC"]+"""
            mkdir -p {output.fastqc} 2> /dev/null && fastqc --quiet --threads {resources.cpus} --outdir {output.fastqc} {input.fastq}/*.fastq.gz &> {log}
        """
###############################################################################
rule fastqscreen_contamination_checking:
    #threads: get_threads('fastqscreen_contamination_checking', 2)
    message: "Fastq-Screen reads contamination checking"
    resources: cpus=1, mem_mb=4000, tim_min=60                  #get_mem('fastqscreen_contamination_checking', 4000)
    params:
        config = CONFIG,
        mapper = MAPPER,
        subset = SUBSET
    input:
        fastq = "resources/reads/"
    output:
        fastqscreen = directory("results/00_Quality_Control/fastq-screen/"),
    log:
        "results/11_Reports/quality/fastq-screen.log"
    shell:
        config["MODULES"]["FASTQSCREEN"]+"\n"+config["MODULES"]["BWA"]+"""
            fastq_screen -q --threads {resources.cpus} --conf {params.config} --aligner {params.mapper} --subset {params.subset} --outdir {output.fastqscreen} {input.fastq}/*.fastq.gz &> {log}
        """

###############################################################################
rule trimmomatic:
    #threads: get_threads('trimmomatic', 8),
    message: "Trimming reads for {wildcards.sample}"
    input:
        r1="resources/reads/{sample}_R1.fastq.gz",
        r2="resources/reads/{sample}_R2.fastq.gz",
        adapters = config["trimmomatic"]["adapters"]["truseq2-pe"]
    output:
        forward_reads   = "results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R1.fastq.gz",
        reverse_reads   = "results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R2.fastq.gz",
        forwardUnpaired = "results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_unpaired_R1.fastq.gz",
        reverseUnpaired = "results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_unpaired_R2.fastq.gz"
    log:
        "results/11_Reports/trimmomatic/{sample}.log"
    params:
        seedMisMatches =            str(config['trimmomatic']['seedMisMatches']),
        palindromeClipTreshold =    str(config['trimmomatic']['palindromeClipTreshold']),
        simpleClipThreshold =      str(config['trimmomatic']['simpleClipThreshold']),
        LeadMinTrimQual =           str(config['trimmomatic']['LeadMinTrimQual']),
        TrailMinTrimQual =          str(config['trimmomatic']['TrailMinTrimQual']),
        windowSize =                str(config['trimmomatic']['windowSize']),
        avgMinQual =                str(config['trimmomatic']['avgMinQual']),
        minReadLen =                str(config['trimmomatic']['minReadLength']),
        phred = 		            str(config["trimmomatic"]["phred"])
    resources: cpus=8, mem_mb=6000, time_min=300                         #get_mem('trimmomatic', 6000)
    shell:
        config["MODULES"]["TRIMMOMATIC"]+"""
            trimmomatic PE -threads {resources.cpus} {params.phred} {input.r1} {input.r2} \ 
            {output.forward_reads} {output.forwardUnpaired} {output.reverse_reads} \
            {output.reverseUnpaired} \
            ILLUMINACLIP:{input.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshold} \
            LEADING:20 \
            TRAILING:3 \
            SLIDINGWINDOW:5:20 \
            AVGQUAL:20 \
            MINLEN:50 \
            &>{log}
        """
###############################################################################
rule trimmed_fastqc:
    #threads: get_threads('fastqc_quality_control', 1)    
    message: "Quality check after trimming"
    resources: cpus=1, mem_mb=4000, time_min=120                            #get_mem('fastqc_quality_control', 4000)
    input:
        forward_reads = expand("results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R1.fastq.gz", sample=SAMPLE),
        reverse_reads = expand("results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R2.fastq.gz", sample=SAMPLE),
    output:
        fastqc = directory("results/00_Quality_Control/trimmed_fastqc/")
    log:
        "results/11_Reports/trimmed_fastqc/trimmed.fastqc.log"
    shell:
        config["MODULES"]["FASTQC"]+"""
            mkdir -p {output.fastqc} 2> /dev/null && \
            fastqc --quiet --threads {resources.cpus} --outdir {output.fastqc} \
            {input.forward_reads} {input.reverse_reads} &> {log}
            """
