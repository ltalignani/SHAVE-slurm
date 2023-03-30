#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 shave.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment pipeline for VEctor 
# Date:                 2022.10.05
# Run:                  snakemake s- workflow/rules/shave.smk --profile slurm/
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
# ENVIRONMENTS #

###############################################################################
# PARAMETERS #
LENGTHc = config["cutadapt"]["length"]              # Cutadapt --minimum-length
TRUSEQ = config["cutadapt"]["kits"]["truseq"]       # Cutadapt --adapter Illumina TruSeq
NEXTERA = config["cutadapt"]["kits"]["nextera"]     # Cutadapt --adapter Illumina Nextera
SMALL = config["cutadapt"]["kits"]["small"]         # Cutadapt --adapter Illumina Small

COMMAND = config["sickle-trim"]["command"]          # Sickle-trim command
ENCODING = config["sickle-trim"]["encoding"]        # Sickle-trim --qual-type
QUALITY = config["sickle-trim"]["quality"]          # Sickle-trim --qual-threshold
LENGTH = config["sickle-trim"]["length"]            # Sickle-trim --length-treshold

# TRIMMOMATIC:
TRUSEQ2_PE: config['trimmomatic']["adapters"]["truseq2-pe"]     # Truseq2-PE adapters fasta
TRUSEQ2_SE: config['trimmomatic']["adapters"]["truseq2-se"]     # Truseq2-SE adapters fasta
TRUSEQ3_PE: config['trimmomatic']["adapters"]["truseq3-pe"]     # Truseq2-PE adapters fasta
TRUSEQ3_PE2: config['trimmomatic']["adapters"]["truseq3-pe-2"]   # Truseq3-PE2 adapters fasta
TRUSEQ3_SE: config['trimmomatic']["adapters"]["truseq3-se"]     # Truseq2-PE adapters fasta

CONFIG = config["fastq-screen"]["config"]           # Fastq-screen --conf
MAPPER = config["fastq-screen"]["aligner"]          # Fastq-screen --aligner
SUBSET = config["fastq-screen"]["subset"]           # Fastq-screen --subset

TRIMMER = config["trimmer"]                         # Trimmers ('sickle' or 'trimmomatic')
ALIGNER = config["aligner"]                         # Aligners ('bwa' or 'bowtie2')
MARKDUP = config["markdup"]			    # Mark Duplicates Program ('picard' or 'samtools')

BWAPATH = config["bwa"]["path"]                     # BWA path to indexes
BT2PATH = config["bowtie2"]["path"]                 # Bowtie2 path to indexes
SENSITIVITY = config["bowtie2"]["sensitivity"]      # Bowtie2 sensitivity preset

REFPATH = config["path"]                            # Path to genomes references
REFERENCE = config["reference"]                     # Genome reference sequence, in fasta format
INDEX = config["index"]                             # Genome reference index, in .fai format
DICTIONARY = config["dictionary"]                  # Genome reference dictionary, made w/ picard CreateSequenceDictionary, in .dict format

ALLELES = config["alleles"]["alleles_target"]          # Alleles against which to genotype (VCF format) 
MINCOV = config["consensus"]["mincov"]              # Minimum coverage, mask lower regions with 'N'
MINAF = config["consensus"]["minaf"]                # Minimum allele frequency allowed
IUPAC = config["consensus"]["iupac"]                # Output variants in the form of IUPAC ambiguity codes

###############################################################################
# FUNCTIONS AND COMMANDS #

############################## O N S T A R T ###################################
onstart:
    shell("mkdir -p Cluster_logs/")

################################## A L L #######################################
rule all:
    input:
        vcf_filtered ="results/05_Variants/merged_filtered/merged_hardfiltered.vcf.gz",

############################### R U L E S #####################################
rule gatk_filter:
    # Aim: Filter variant calls based on INFO and/or FORMAT annotations.
    # Use: gatk VariantFiltration \
    # -R reference.fasta \
    # -V input.vcf.gz \
    # -O output.vcf.gz
    # --filter-name "my_filter1" \
    # --filter-expression "AB < 0.2" \
    # --filter-name "my_filter2" \
    # --filter-expression "MQ0 > 50"
    message:
        "VariantFiltration Hard-filtering"
    input:
        ref=REFPATH+REFERENCE,  #"resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
        vcf="results/05_Variants/merged_raw/merged.vcf.gz",
    output:
        vcf="results/05_Variants/merged_filtered/merged_hardfiltered.vcf.gz",
    params:
        filters={"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"},
        extra="",
        java_opts="",
    resources: cpus=8, mem_mb=12000, time_min=120
    log:
        "results/11_Reports/variantfiltration/merged_hardfiltered.log",
    shell:
        config["MODULES"]["GATK4"]+"""
            gatk VariantFiltration -R {input.ref} -V {input.vcf} -O {output.vcf} --filterExpression {params.filters} --filterName 'my_filters'
        """
