#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 shave_polish_3.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment pipeline for VEctor 
# Date:                 2022.10.05
# Run:                  snakemake -s workflow/rules/shave_polish_3.smk --profile slurm/
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
        qualimap = expand("results/00_Quality_Control/qualimap/{sample}_bwa/qualimapReport.html", sample=SAMPLE),
        index_post_realign = expand("results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bai", sample=SAMPLE),
        fixmateinformation = expand("results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam", sample=SAMPLE),

############################### R U L E S #####################################
rule fixmateinformation:
    # Aim: This tool ensures that all mate-pair information is in sync between each read and its mate pair.
    #      If no #OUTPUT file is supplied then the output is written to a temporary file and then copied over
    #      the #INPUT file (with the original placed in a .old file.)
    # Use: picard.jar FixMateInformation \
    #      -I input.bam \
    #      -O fixed_mate.bam \
    #      --ADD_MATE_CIGAR true
    message:
        "Picard FixMateInformation"
    input:
        realigned = "results/04_Polishing/realigned/{sample}_bwa_md_realigned.bam",
    output:
        fixed = "results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam",
    resources: cpus=8, mem_mb=4000, time_min=120
    log:
        "results/11_Reports/fixmateinformation/{sample}_bwa_realigned_fixed.log",
    shell:
        config["MODULES"]["PICARDTOOLS"]+"""
            picard FixMateInformation -I {input.realigned} -O {output.fixed} --ADD_MATE_CIGAR true &> {log}
        """

###############################################################################
rule samtools_index_post_realign:
    # Aim: indexing marked as duplicate BAM file
    # Use: samtools index -@ [THREADS] -b [MARKDUP.bam] [INDEX.bai]
    message:
        "SamTools indexing realigned fixed-mate sorted BAM file {wildcards.sample} sample for Picard ValidateSamFile"
    resources: cpus=6, mem_mb=4000, time_min=120
    input:
        fixedbam = rules.fixmateinformation.output.fixed                                #"results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam",
    output:
        index = "results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bai",
    log:
        "results/11_Reports/samtools/{sample}_bwa_realigned_fixed_indexed.log",
    shell:
        config["MODULES"]["SAMTOOLS"]+"""
            samtools index -@ {threads} -b {input.fixedbam} {output.index} &> {log}
        """

###############################################################################
rule qualimap:
    # Aim: Qualimap is a platform-independent application written in Java and R that provides both a Graphical User Interface (GUI) and a 
    # command-line interface to facilitate the quality control of alignment sequencing data. Shortly, Qualimap:
    #       1. Examines sequencing alignment data according to the features of the mapped reads and their genomic properties
    #       2. Povides an overall view of the data that helps to to the detect biases in the sequencing and/or mapping of the data and eases decision-making for further analysis.
    #
    # Use: qualimap bamqc -bam {input.bam} \ 
    #       -c \                        Paint chromosome limits inside charts
    #       -nt {threads} \             
    #       -outdir {output.report} \   Output directory for HTML report (default value is report.html)
    #       -outformat PDF \            Format of the ouput report (PDF or HTML, default is HTML)
    #       -sd                         Activate this option to skip duplicate alignments from the analysis                   
    input:
        bam = rules.fixmateinformation.output.fixed,                                 #"results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam",
        index = rules.samtools_index_post_realign.output.index                      # not used in the command, but it's here so snakemake knows to run the rule after the indexing
    output:
        protected("results/00_Quality_Control/qualimap/{sample}_bwa/qualimapReport.html"),
        protected("results/00_Quality_Control/qualimap/{sample}_bwa/raw_data_qualimapReport/genome_fraction_coverage.txt"),
        protected("results/00_Quality_Control/qualimap/{sample}_bwa/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt"),
        protected("results/00_Quality_Control/qualimap/{sample}_bwa/genome_results.txt"),
        protected("results/00_Quality_Control/qualimap/{sample}_bwa/raw_data_qualimapReport/coverage_histogram.txt")
    params:
        outdir = "results/00_Quality_Control/qualimap/{sample}/"
    resources: cpus=6, mem_mb=8000, time_min=120
    log:
        "results/11_Reports/qualimap/logs/{sample}_bwa_qualimap.log" 
    shell:
        config["MODULES"]["QUALIMAP"]+"""
            unset DISPLAY && qualimap bamqc -bam {input.bam} -nt {threads} --java-mem-size=8G -outdir {params.outdir} 
        """
#qualimap bamqc -bam {input.bam} -c -nt {threads} --java-mem-size={resources.mem_mb}M -outdir results/00_Quality_Control/qualimap/{wildcards.sample} -sd &> {log}