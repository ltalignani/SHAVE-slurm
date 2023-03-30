#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 shave_polishing.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment pipeline for VEctor 
# Date:                 2022.10.05
# Run:                  snakemake -s workflow/rules/shave_polishing.smk --profile slurm/
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
#     if "__default__" in cluster_config and "memem-per-cpum" in cluster_config["__default__"]:
#         return int(cluster_config["__default__"]["mem-per-cpu"])
#     return default

############################## O N S T A R T ###################################
onstart:
    shell("mkdir -p Cluster_logs/")

################################## A L L #######################################
rule all:
    input:
        fix                 = expand("results/02_Mapping/{sample}_bwa_sorted-mark-dup-fx.bam", sample=SAMPLE),
        index               = expand("results/02_Mapping/{sample}_bwa_sorted-mark-dup-fx.bai", sample=SAMPLE),
        intervals           = expand("results/04_Polishing/{sample}_bwa.intervals", sample=SAMPLE),
        bed                 = expand("results/04_Polishing/{sample}_bwa_realignertargetcreator.bed", sample=SAMPLE),
        bam                 = expand("results/04_Polishing/realigned/{sample}_bwa_md_realigned.bam", sample=SAMPLE),
        bai                 = expand("results/04_Polishing/realigned/{sample}_bwa_md_realigned.bai", sample=SAMPLE),
        qualimap            = expand("results/00_Quality_Control/qualimap/{sample}_bwa/qualimapReport.html", sample=SAMPLE),
        index_post_realign  = expand("results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bai", sample=SAMPLE),
        fixmateinformation  = expand("results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam", sample=SAMPLE),
        multiqc             = "multiqc_report.html",
        quickcheck          = "bad_bams.fofn",
        check               = expand("results/00_Quality_Control/validatesamfile/{sample}_bwa_md_realigned_fixed_ValidateSam.txt", sample=SAMPLE),
        flagstat            = expand("results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed_bam.flagstat.txt", sample=SAMPLE),
        idxstats            = expand("results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed.idxstats.txt", sample=SAMPLE),
        stats               = expand("results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed_stats.txt", sample=SAMPLE),     


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
    envmodules:
        "picard/2.23.5"
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

################################################################################
rule indelrealigner:
    # Aim:  Mappers cannot “see” indels near ends of reads because mismatches are “cheaper” than a gap in this context.
    #       IndelRealigner takes a coordinate-sorted and indexed BAM and a target intervals file generated by RealignerTargetCreator.
    #       IndelRealigner then performs local realignment on reads coincident with the target intervals using consenses
    #       from indels present in the original alignment.
    # Use:  java -Xmx16G -Djava.io.tmpdir=/tmp -jar GenomeAnalysisTK.jar \
    #       -T IndelRealigner \
    #       -R human_g1k_v37_decoy.fasta \
    #       -targetIntervals realignertargetcreator.intervals \
    #       -known INDEL_chr10_1Mb_b37_1000G_phase3_v4_20130502.vcf \
    #       -I 7156_snippet.bam \
    #       -o 7156_snippet_indelrealigner.bam
    message:
        "Indel realignment"
    input:
        bam="results/02_Mapping/{sample}_bwa_sorted-mark-dup-fx.bam",
        bai="results/02_Mapping/{sample}_bwa_sorted-mark-dup-fx.bai",
        ref=REFPATH+REFERENCE,                                                            #"resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
        fai="resources/genomes/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa.fai",
        dict="resources/genomes/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.dict",
        target_intervals="results/04_Polishing/{sample}_bwa.intervals"
    output:
        bam= "results/04_Polishing/realigned/{sample}_bwa_md_realigned.bam",
        bai= "results/04_Polishing/realigned/{sample}_bwa_md_realigned.bai",
    benchmark:
        "benchmarks/indelrealigner/{sample}_bwa.tsv",
    log:
        "results/11_Reports/indelrealigner/{sample}_bwa.log",
    params:
        extra="--defaultBaseQualities 20 --filter_reads_with_N_cigar",  # optional
    resources: cpus=12, mem_mb=12000, time_min=120
    shell:
        config["MODULES"]["GATK3"]+"""
            gatk3 -T IndelRealigner -I {input.bam} -R {input.ref} --targetIntervals {input.target_intervals} {params.extra} --out {output.bam}
        """
        
        #"java -jar GenomeAnalysisTK.jar -T IndelRealigner -I {input.bam} -R {input.ref} --targetIntervals {input.target_intervals} {params.extra} --out {output.bam}"

################################################################################
rule awk_intervals_for_IGV:
    # Aim: View intervals on IGV
    # Use: awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' \
    #      realignertargetcreator.intervals > realignertargetcreator.bed
    message:
        "Awk IGV intervals visualization for {wildcards.sample} sample "
    input:
        intervals="results/04_Polishing/{sample}_bwa.intervals",
    resources: cpus=1, mem_mb=4000, time_min=120
    params:
        cmd = r"""'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}'"""
    output:
        bed = "results/04_Polishing/{sample}_bwa_realignertargetcreator.bed",
    log:
        "results/11_Reports/awk/{sample}_bwa_intervals_for_IGV.log"
    shell:
        "awk -F '[:-]' {params.cmd} {input.intervals} 1> {output.bed} 2> {log}"

###############################################################################
rule realignertargetcreator:
    # Aim:      RealignerTargetCreator identify what regions need to be realigned.
    #           Local realignment around indels. Takes a coordinate-sorted and indexed BAM and a VCF of known indels and creates a target intervals file.
    # Use:      gatk3 -T RealignerTargetCreator \
    #           -R human_g1k_v37_decoy.fasta \
    #           -L 10:96000000-97000000 \
    #           -known INDEL_chr10_1Mb_b37_1000G_phase3_v4_20130502.vcf \
    #           -I 7156_snippet.bam \
    #           --out 7156_realignertargetcreator.intervals
    message:
        "RealignerTargetCreator creates a target intervals file for indel realignment"
    input:
        bam="results/02_Mapping/{sample}_bwa_sorted-mark-dup-fx.bam",
        bai="results/02_Mapping/{sample}_bwa_sorted-mark-dup-fx.bai",
        ref=REFPATH+REFERENCE,
        fai="resources/genomes/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa.fai",
        dict="resources/genomes/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.dict",
    output:
        "results/04_Polishing/{sample}_bwa.intervals",
    benchmark:
        "benchmarks/realignertargetcreator/{sample}_bwa.tsv",
    log:
        "results/11_Reports/realignertargetcreator/{sample}_bwa.log",
    params:
        extra="--defaultBaseQualities 20 --filter_reads_with_N_cigar",  # optional
    resources: cpus=8, mem_mb=12000, time_min=120
    shell:
        config["MODULES"]["GATK3"]+"""
            gatk3 -T RealignerTargetCreator --num_threads {threads} -R {input.ref} -I {input.bam} {params.extra} --out {output}
        """
        
        #"java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator --num_threads {threads} -R {input.ref} -I {input.bam} {params.extra} --out {output} "

###############################################################################
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
        
###############################################################################
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
    output:
        "multiqc_report.html"
    resources: cpus=4, mem_mb=8000, time_min=120
    log:
        "results/11_Reports/multiqc/multiqc.log"
    shell:
        config["MODULES"]["MULTIQC"]+"""
        multiqc {input} -o results/00_Quality_Control/MULTIQC/ -n {output} > {log} 2>&1 
        """
