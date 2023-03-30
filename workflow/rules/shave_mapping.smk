#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 shave_mapping.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment pipeline for VEctor 
# Date:                 2022.10.05
# Run:                  snakemake -s workflow/rules/shave_mapping.smk --profile slurm/
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

############################## O N S T A R T ###################################
onstart:
    shell("mkdir -p Cluster_logs/")

################################## A L L #######################################
rule all:
    input:
        mapped = expand("results/02_Mapping/{sample}_bwa-mapped.sam", sample= SAMPLE),

############################### R U L E S #####################################
rule bwa_mapping:
    # Aim: reads mapping against reference sequence
    # Use: bwa mem -t [THREADS] -x [REFERENCE] [FWD_R1.fq] [REV_R2.fq] 1> [MAPPED.sam]
    message:
        "BWA-MEM mapping {wildcards.sample} sample reads against reference genome sequence"
    resources: cpus=16, mem_mb=16000, time_min=600
    params:
        ref = REFPATH+REFERENCE,
        extra = r"'@RG\tID:{sample}\tSM:{sample}\tCN:SC\tPL:ILLUMINA'" # Manage ReadGroup
    input:
        fwdreads = "results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R1.fastq.gz",
        revreads = "results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R2.fastq.gz"
    output:
        mapped = "results/02_Mapping/{sample}_bwa-mapped.sam",
    benchmark:
        "benchmarks/bwa/{sample}.tsv"
    log:
        "results/11_Reports/bwa/{sample}.log"
    shell:
        config["MODULES"]["BWA"]+"""
            bwa mem -M -T 0 -t {threads} -v 1 -R {params.extra} {params.ref} {input.fwdreads} {input.revreads} 1> {output.mapped} 2> {log}
        """