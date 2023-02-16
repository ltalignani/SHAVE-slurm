#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 shave_unifiedgenotyper.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment pipeline for VEctor 
# Date:                 2022.10.05
# Run:                  snakemake -s workflow/rules/shave_unifiedgenotyper.smk \
#                       --cluster-config cluster.yaml --configfile config.yaml 
#                       --rerun-incomplete --cores 24 --latency-wait 600
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

################################## A L L #######################################
rule all:
    input:
        table = "results/05_Variants/merged_raw/merged.table",
        indexvcf = expand("results/05_Variants/{sample}_bwa.vcf.gz.tbi", sample=SAMPLE),
        combinegvcfs = "results/05_Variants/merged_raw/merged.vcf.gz",
        bgzip_vcfs = expand("results/05_Variants/{sample}_bwa.vcf.gz", sample=SAMPLE),
        vcf = expand("results/05_Variants/{sample}_bwa.vcf", sample=SAMPLE),

################################ R U L E S #####################################
rule unifiedgenotyper:
    # Aim:  Call variants in sequence data. The following parameters comes from the MalariaGEN
    # Use:  "gatk3 -T UnifiedGenotyper "                    # Genome Analysis Tool Kit - Broad Institute UnifiedGenotyper
    #    "-nct {threads} "                               # -nct / --num_cpu_threads_per_data_thread controls the number of CPU threads allocated to each data thread
    #    "-I {input.bam} "                               # Input indel realigned BAM file
    #    "-R {input.ref} "                               # Reference sequence in fasta format
    #    "--out {output.vcf} "                           # Output VCF
    #    "--genotype_likelihoods_model BOTH "            # Genotype likelihoods calculation model to employ -- BOTH is the default option, while INDEL is also available for calling indels and SNP is available for calling SNPs only (SNP|INDEL|BOTH)
    #    "--genotyping_mode DISCOVERY "                  # Should we output confident genotypes (i.e. including ref calls) or just the variants? (DISCOVERY|GENOTYPE_GIVEN_ALLELES)
    #    "--heterozygosity 0.015 "                       # Heterozygosity value used to compute prior likelihoods for any locus
    #    "--heterozygosity_stdev 0.05 "                  # Standard deviation of heterozygosity for SNP and indel calling
    #    "--indel_heterozygosity 0.001 "                 # Heterozygosity for indel calling
    #    "--downsampling_type BY_SAMPLE "                # Type of reads downsampling to employ at a given locus. Reads will be selected randomly to be removed from thepile based on the method described here (NONE|ALL_READS| BY_SAMPLE) given locus
    #    "-dcov 250 "                                    # downsampling coverage
    #    "--output_mode EMIT_ALL_SITES "                 # Should we output confident genotypes (i.e. including ref calls) or just the variants? (EMIT_VARIANTS_ONLY|EMIT_ALL_CONFIDENT_SITES|EMIT_ALL_SITES)
    #    "--min_base_quality_score 17 "                  # Minimum base quality required to consider a base for calling
    #    "-stand_call_conf 0.0 "                         # standard min confidence-threshold for calling
    #    "-contamination 0.0 "                           # Define the fraction of contamination in sequence data (for all samples) to aggressively remove.
    #    "-A DepthPerAlleleBySample "                    #
    #    "-A RMSMappingQuality "                         # MQ: needed as decision tools for hard-filtering. Compares the mapping qualities of the reads supporting the reference allele and the alternate allele. positive value means the mapping qualities of the reads supporting the alternate allele are higher than those supporting the reference allele
    #    "-A Coverage "                                  #
    #    "-A FisherStrand "                              # FS: needed as decision tools for hard-filtering. Strand Bias tells us whether the alternate allele was seen more or less often on the forward or reverse strand than the reference allele. measure strand bias (a type of sequencing bias in which one DNA strand is favored over the other, which can result in incorrect evaluation of the amount of evidence observed for one allele vs. the other
    #    "-A StrandOddsRatio "                           # SOR: needed as decision tools for hard-filtering. created because FS tends to penalize variants that occur at the ends of exons. Reads at the ends of exons tend to only be covered by reads in one direction
    #    "-A BaseQualityRankSumTest "                    #
    #    "-A MappingQualityRankSumTest "                 # MQRankSum: needed as decision tools for hard-filtering
    #    "-A QualByDepth "                               # QD: needed as decision tools for hard-filtering. Intended to normalize the variant quality in order to avoid inflation caused when there is deep coverage. This is the variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples. better to use QD than either QUAL or DP directly.
    #    "-A ReadPosRankSumTest "                        # ReadPosRankSum: needed as decision tools for hard-filtering. z-approximation from the Rank Sum Test for site position within reads. Compares whether the positions of the reference and alternate alleles are different within the reads.
    #    "-XA ExcessHet "                                #
    #    "-XA InbreedingCoeff "                          #
    #    "-XA MappingQualityZero "                       #
    #    "-XA HaplotypeScore "                           #
    #    "-XA SpanningDeletions "                        #
    #    "-XA ChromosomeCounts "                         #
    #    " > {log} 2>&1"
    #threads: 24             #get_threads('unifiedgenotyper', 24) <- not necessary with simple slurm profile
    message:
        "UnifiedGenotyper calling SNVs"
    input:
        bam = "results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam",
        ref = "resources/genomes/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa",
        index = "results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bai",
        #alleles = ALLELES
    output:
        vcf="results/05_Variants/{sample}_bwa.vcf"
    log:
        "results/11_Reports/unifiedgenotyper/{sample}_bwa.log"
    benchmark:
        "benchmarks/unifiedgenotyper/{sample}_bwa.tsv"
    resources: cpus=24, mem_mb=16000, time_min=1200       #get_mem('unifiedgenotyper', 12000)  <- not necessary with simple slurm profile
    params: partition = 'long'
    shell:
        config["MODULES"]["GATK3"]+"""
            gatk3 -T UnifiedGenotyper -nct {resources.cpus} -I {input.bam} -R {input.ref} --out {output.vcf} --genotype_likelihoods_model BOTH \
            --genotyping_mode DISCOVERY --heterozygosity 0.015 --heterozygosity_stdev 0.05 --indel_heterozygosity 0.001 --downsampling_type BY_SAMPLE \
            -dcov 250 --output_mode EMIT_ALL_SITES --min_base_quality_score 17 -stand_call_conf 0.0 -contamination 0.0 -A DepthPerAlleleBySample \
            -A RMSMappingQuality -A Coverage -A FisherStrand -A StrandOddsRatio -A BaseQualityRankSumTest -A MappingQualityRankSumTest -A QualByDepth \
            -A ReadPosRankSumTest -XA ExcessHet -XA InbreedingCoeff -XA MappingQualityZero -XA HaplotypeScore -XA SpanningDeletions -XA ChromosomeCounts \
            &> {log}
        """

###############################################################################
rule bgzip_vcfs:
    #threads: 6                      #get_threads('bgzip_vcfs', 6),
    input:
        rules.unifiedgenotyper.output.vcf       #"results/05_Variants/{sample}_bwa.vcf",
    output:
        gzip = "results/05_Variants/{sample}_bwa.vcf.gz",
    resources: cpus=1, mem_mb=4000, time_min=300               #get_mem('bgzip_vcfs', 4000),
    params: partition = 'fast'
    log:
        "results/11_Reports/bgzip/{sample}_bwa.vcf.gz.log",
    shell:
        config["MODULES"]["HTSLIB"]+"""
            bgzip -c --threads {threads} {input} > {output.gzip} {log}
        """

###############################################################################
rule indexfeaturefile:
    #threads: 4    
    message:
        "Indexing VCFs file for CombineGVCFs"
    input:
        vcf = rules.bgzip_vcfs.output.gzip                  #"results/05_Variants/{sample}_bwa.vcf.gz",
    output:
        indexvcf = "results/05_Variants/{sample}_bwa.vcf.gz.tbi",
    resources: cpus=1, mem_mb=4000, time_min=120
    params: partition = 'fast'
    log:
        "results/11_Reports/indexfeaturefile/{sample}_bwa.indexvcf.tbi.log",
    shell:
        config["MODULES"]["GATK4"]+"""
            gatk IndexFeatureFile -I {input.vcf} -O {output.indexvcf} > {log} 2>&1 || true
        """

###############################################################################
rule combinegvcfs:
    #threads: 8    
    message:
        "GATK CombineGVCFs merge VCFs"
    input:
        gvcfs= expand("results/05_Variants/{sample}_bwa.vcf.gz", sample=SAMPLE),
        ref=REFPATH+REFERENCE,
    output:
        gvcf="results/05_Variants/merged_raw/merged.vcf.gz",
    log:
        "results/11_Reports/combinegvcfs/combinegvcfs.log",
    resources: cpus=1, mem_mb=8000, time_min=120           #get_mem('combinegvcfs', 12000)
    params: partition = 'fast'
    shell:
        config["MODULES"]["GATK4"]+"""
            gatk CombineGVCFs -R {input.ref} -I {input.gvcfs} -O {output.gvcf} &> {log} 2>&1 || true
        """

###############################################################################
rule variantstotable:
    # Aim:  Call variants in sequence data. The following parameters comes from the MalariaGEN
    # Use:  gatk VariantsToTable\
    #       -V input.vcf \
    #       -F CHROM -F POS -F TYPE -GF AD \
    #       -O output.table
    #threads: 2
    message:
        "VariantsToTable"
    input:
        vcf = rules.combinegvcfs.output.gvcf            #"results/05_Variants/merged_raw/merged.vcf.gz",
    output:
        table = "results/05_Variants/merged_raw/merged.table",
    resources: cpus=1, mem_mb=4000, tim_min=60
    params: partition = 'fast'
    log:
        "results/11_Reports/variantstotable/merged_table.log",
    shell:
        config["MODULES"]["GATK4"]+"""
            gatk VariantsToTable -V {input.vcf} -O {output.table} -F CHROM -F POS -F QUAL -F NS -F DP -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum -GF AD 2> {log}
        """
