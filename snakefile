#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 snakefile
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment pipeline for VEctor 
# Date:                 2022.10.05
# Run:                  snakemake --profile slurm/
# Latest modification:  2023.01.25
# Done:                 Added cluster-config, get_threads(),

###############################################################################
# CONFIGURATION FILES #
configfile: "config/config.yaml"
cluster_config: "cluster.yaml"

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
DICTIONARY = config["dictionary"]                   # Genome reference dictionary, made w/ picard CreateSequenceDictionary, in .dict format

ALLELES = config["alleles"]["alleles_target"]          # Alleles against which to genotype (VCF format) 
MINCOV = config["consensus"]["mincov"]              # Minimum coverage, mask lower regions with 'N'
MINAF = config["consensus"]["minaf"]                # Minimum allele frequency allowed
IUPAC = config["consensus"]["iupac"]                # Output variants in the form of IUPAC ambiguity codes

###############################################################################
# FUNCTIONS AND COMMANDS #

################################ O N S T A R T #################################
onstart:
    shell("mkdir -p Cluster_logs/")

###################### R U L E   D E C L A R A T I O N #########################

include_prefix=os.getcwd()+"/"

include:
    include_prefix + "workflow/rules/shave_trimmomatic.smk"
include:
    include_prefix + "workflow/rules/shave_mapping.smk"
include:
    include_prefix + "workflow/rules/shave_markdup.smk"
include:
    include_prefix + "workflow/rules/shave_polish1.smk"
include:
    include_prefix + "workflow/rules/shave_polish2.smk"
include:
    include_prefix + "workflow/rules/shave_polish3.smk"
include:
    include_prefix + "workflow/rules/shave_polish4.smk"
include:
    include_prefix + "workflow/rules/shave_unifiedgenotyper.smk"
include:
    include_prefix + "workflow/rules/shave_vcf_hf.smk"

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
        mapped          =   expand("results/02_Mapping/{sample}_bwa-mapped.sam", sample= SAMPLE),
        bam             =   expand("results/02_Mapping/{sample}_bwa_sorted-mark-dup.bam", sample=SAMPLE),
        metrics         =   expand("results/02_Mapping/{sample}_bwa_sorted-mark-dup_metrics.txt", sample=SAMPLE),    
        fix             =   expand("results/02_Mapping/{sample}_bwa_sorted-mark-dup-fx.bam", sample=SAMPLE),
        index           =   expand("results/02_Mapping/{sample}_bwa_sorted-mark-dup-fx.bai", sample=SAMPLE),
        intervals       =   expand("results/04_Polishing/{sample}_bwa.intervals", sample=SAMPLE),
        bed             =   expand("results/04_Polishing/{sample}_bwa_realignertargetcreator.bed", sample=SAMPLE),
        bam             =   expand("results/04_Polishing/realigned/{sample}_bwa_md_realigned.bam", sample=SAMPLE),
        bai             =   expand("results/04_Polishing/realigned/{sample}_bwa_md_realigned.bai", sample=SAMPLE),
        qualimap        =   expand("results/00_Quality_Control/qualimap/{sample}_bwa/qualimapReport.html", sample=SAMPLE),
        index_post_realign  = expand("results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bai", sample=SAMPLE),
        fixmateinformation  = expand("results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam", sample=SAMPLE),
        multiqc             = "multiqc_report.html",
        check               = expand("results/00_Quality_Control/validatesamfile/{sample}_bwa_md_realigned_fixed_ValidateSam.txt", sample=SAMPLE),
        flagstat            = expand("results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed_bam.flagstat.txt", sample=SAMPLE),
        idxstats            = expand("results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed.idxstats.txt", sample=SAMPLE),
        stats               = expand("results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed_stats.txt", sample=SAMPLE),
        vcf_filtered        = "results/05_Variants/merged_filtered/merged_hardfiltered.vcf.gz",                      