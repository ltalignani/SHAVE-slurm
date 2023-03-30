#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 shave_unifiedgenotyper.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment pipeline for VEctor 
# Date:                 2022.10.05
# Run:                  snakemake -s workflow/rules/shave_unifiedgenotyper.smk \
#                       --cluster-config cluster.yaml --configfile config.yaml 
#                       --rerun-incomplete --cores 24 --latency-wait 1200
# Latest modification:  2023.03.13
# Done:                 Added VCF stats, VCF reports, send all bams files in UG

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

################################ O N S T A R T #################################
onstart:
    shell("ls results/04_Polishing/realigned/*_fixed.bam > bam.list")

################################## A L L #######################################
rule all:
    input:
        report = "results/report_vcf_X.html",
        report2L = "results/report_vcf_2L.html",
        report2R = "results/report_vcf_2R.html",
        report3L= "results/report_vcf_3L.html",
        report3R = "results/report_vcf_3R.html",
        table = "results/05_Variants/vcf_X.table",
        table2L = "results/05_Variants/vcf_2L.table",
        table2R = "results/05_Variants/vcf_2R.table",
        table3L = "results/05_Variants/vcf_3L.table",
        table3R = "results/05_Variants/vcf_3R.table",
        indexvcf = "results/05_Variants/variants_X.vcf.gz.tbi",
        indexvcf_2L = "results/05_Variants/variants_2L.vcf.gz.tbi",
        indexvcf_2R = "results/05_Variants/variants_2R.vcf.gz.tbi",
        indexvcf_3L = "results/05_Variants/variants_3L.vcf.gz.tbi",
        indexvcf_3R = "results/05_Variants/variants_3R.vcf.gz.tbi",
        bgzip_vcfs = expand("results/05_Variants/variants_X.vcf.gz", sample=SAMPLE),
        bgzip_vcfs2L = expand("results/05_Variants/variants_2L.vcf.gz", sample=SAMPLE),
        bgzip_vcfs2R = expand("results/05_Variants/variants_2R.vcf.gz", sample=SAMPLE),
        bgzip_vcfs3L = expand("results/05_Variants/variants_3L.vcf.gz", sample=SAMPLE),
        bgzip_vcfs3R = expand("results/05_Variants/variants_3R.vcf.gz", sample=SAMPLE),
        vcf = expand("results/05_Variants/variants_X.vcf", sample=SAMPLE),
        vcf2L = expand("results/05_Variants/variants_2L.vcf", sample=SAMPLE),
        vcf2R = expand("results/05_Variants/variants_2R.vcf", sample=SAMPLE),
        vcf3L = expand("results/05_Variants/variants_3L.vcf", sample=SAMPLE),
        vcf3R = expand("results/05_Variants/variants_3R.vcf", sample=SAMPLE),

################################ R U L E S #####################################
rule unifiedgenotyper:
    # Aim:  Call variants in sequence data. The following parameters comes from the MalariaGEN
    # Use:  "gatk3 -T UnifiedGenotyper "                 # Genome Analysis Tool Kit - Broad Institute UnifiedGenotyper
    #    "--num_threads {threads} "                      # -nt / --num_threads 
    #    "-nct {threads} "                               # -nct / --num_cpu_threads_per_data_thread controls the number of CPU threads allocated to each data thread
    #    "-I {input.bam} "                               # Input indel realigned BAM file. Possibility to add a .list file containing path to each bam file (one path per line)
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
    #    " > {log} &>&1"
    message:
        "UnifiedGenotyper calling SNVs for chromosome X"
    input:
        bam = "bam.list",
        ref = "resources/genomes/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa",
    output:
        vcf="results/05_Variants/variants_X.vcf"
    log:
        "results/11_Reports/unifiedgenotyper/variants_X.log"
    benchmark:
        "benchmarks/unifiedgenotyper/variants_ug.tsv"
    resources: cpus=1, mem_mb=112000, time_min=14400       
    params: partition = 'long'
    shell:
        config["MODULES"]["GATK3"]+"""
            gatk3 -T UnifiedGenotyper --num_threads {resources.cpus} --num_cpu_threads_per_data_thread {resources.cpus} -I {input.bam} -R {input.ref} \
            -L X --out {output.vcf} --genotype_likelihoods_model BOTH \
            --genotyping_mode DISCOVERY --heterozygosity 0.015 --heterozygosity_stdev 0.05 --indel_heterozygosity 0.001 --downsampling_type BY_SAMPLE \
            -dcov 250 --output_mode EMIT_ALL_SITES --min_base_quality_score 17 -stand_call_conf 0.0 -contamination 0.0 -A DepthPerAlleleBySample \
            -A RMSMappingQuality -A Coverage -A FisherStrand -A StrandOddsRatio -A BaseQualityRankSumTest -A MappingQualityRankSumTest -A QualByDepth \
            -A ReadPosRankSumTest -XA ExcessHet -XA InbreedingCoeff -XA MappingQualityZero -XA HaplotypeScore -XA SpanningDeletions -XA ChromosomeCounts \
            &> {log}
        """

###############################################################################
rule unifiedgenotyper_2L:
    message:
        "UnifiedGenotyper calling SNVs for chromosome 2L"
    input:
        bam = "bam.list",
        ref = "resources/genomes/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa",
    output:
        vcf="results/05_Variants/variants_2L.vcf"
    log:
        "results/11_Reports/unifiedgenotyper/variants_2L.log"
    benchmark:
        "benchmarks/unifiedgenotyper/variants_2L.tsv"
    resources: cpus=1, mem_mb=112000, time_min=14400       
    params: partition = 'long'
    shell:
        config["MODULES"]["GATK3"]+"""
            gatk3 -T UnifiedGenotyper --num_threads {resources.cpus} --num_cpu_threads_per_data_thread {resources.cpus} -I {input.bam} -R {input.ref} \
            -L 2L --out {output.vcf} --genotype_likelihoods_model BOTH \
            --genotyping_mode DISCOVERY --heterozygosity 0.015 --heterozygosity_stdev 0.05 --indel_heterozygosity 0.001 --downsampling_type BY_SAMPLE \
            -dcov 250 --output_mode EMIT_ALL_SITES --min_base_quality_score 17 -stand_call_conf 0.0 -contamination 0.0 -A DepthPerAlleleBySample \
            -A RMSMappingQuality -A Coverage -A FisherStrand -A StrandOddsRatio -A BaseQualityRankSumTest -A MappingQualityRankSumTest -A QualByDepth \
            -A ReadPosRankSumTest -XA ExcessHet -XA InbreedingCoeff -XA MappingQualityZero -XA HaplotypeScore -XA SpanningDeletions -XA ChromosomeCounts \
            &> {log}
        """

###############################################################################
rule unifiedgenotyper_2R:
    message:
        "UnifiedGenotyper calling SNVs for chromosome 2R"
    input:
        bam = "bam.list",
        ref = "resources/genomes/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa",
    output:
        vcf="results/05_Variants/variants_2R.vcf"
    log:
        "results/11_Reports/unifiedgenotyper/variants_2R.log"
    benchmark:
        "benchmarks/unifiedgenotyper/variants_2R.tsv"
    resources: cpus=1, mem_mb=112000, time_min=14400       
    params: partition = 'long'
    shell:
        config["MODULES"]["GATK3"]+"""
            gatk3 -T UnifiedGenotyper --num_threads {resources.cpus} --num_cpu_threads_per_data_thread {resources.cpus} -I {input.bam} -R {input.ref} \
            -L 2R --out {output.vcf} --genotype_likelihoods_model BOTH \
            --genotyping_mode DISCOVERY --heterozygosity 0.015 --heterozygosity_stdev 0.05 --indel_heterozygosity 0.001 --downsampling_type BY_SAMPLE \
            -dcov 250 --output_mode EMIT_ALL_SITES --min_base_quality_score 17 -stand_call_conf 0.0 -contamination 0.0 -A DepthPerAlleleBySample \
            -A RMSMappingQuality -A Coverage -A FisherStrand -A StrandOddsRatio -A BaseQualityRankSumTest -A MappingQualityRankSumTest -A QualByDepth \
            -A ReadPosRankSumTest -XA ExcessHet -XA InbreedingCoeff -XA MappingQualityZero -XA HaplotypeScore -XA SpanningDeletions -XA ChromosomeCounts \
            &> {log}
        """

###############################################################################
rule unifiedgenotyper_3L:
    message:
        "UnifiedGenotyper calling SNVs for chromosome 3L"
    input:
        bam = "bam.list",
        ref = "resources/genomes/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa",
    output:
        vcf="results/05_Variants/variants_3L.vcf"
    log:
        "results/11_Reports/unifiedgenotyper/variants_3L.log"
    benchmark:
        "benchmarks/unifiedgenotyper/variants_3L.tsv"
    resources: cpus=1, mem_mb=112000, time_min=14400       
    params: partition = 'long'
    shell:
        config["MODULES"]["GATK3"]+"""
            gatk3 -T UnifiedGenotyper --num_threads {resources.cpus} --num_cpu_threads_per_data_thread {resources.cpus} -I {input.bam} -R {input.ref} \
            -L 3L --out {output.vcf} --genotype_likelihoods_model BOTH \
            --genotyping_mode DISCOVERY --heterozygosity 0.015 --heterozygosity_stdev 0.05 --indel_heterozygosity 0.001 --downsampling_type BY_SAMPLE \
            -dcov 250 --output_mode EMIT_ALL_SITES --min_base_quality_score 17 -stand_call_conf 0.0 -contamination 0.0 -A DepthPerAlleleBySample \
            -A RMSMappingQuality -A Coverage -A FisherStrand -A StrandOddsRatio -A BaseQualityRankSumTest -A MappingQualityRankSumTest -A QualByDepth \
            -A ReadPosRankSumTest -XA ExcessHet -XA InbreedingCoeff -XA MappingQualityZero -XA HaplotypeScore -XA SpanningDeletions -XA ChromosomeCounts \
            &> {log}
        """

###############################################################################
rule unifiedgenotyper_3R:
    message:
        "UnifiedGenotyper calling SNVs for chromosome 3R"
    input:
        bam = "bam.list",
        ref = "resources/genomes/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa",
    output:
        vcf="results/05_Variants/variants_3R.vcf"
    log:
        "results/11_Reports/unifiedgenotyper/variants_3R.log"
    benchmark:
        "benchmarks/unifiedgenotyper/variants_3R.tsv"
    resources: cpus=1, mem_mb=112000, time_min=14400       
    params: partition = 'long'
    shell:
        config["MODULES"]["GATK3"]+"""
            gatk3 -T UnifiedGenotyper --num_threads {resources.cpus} --num_cpu_threads_per_data_thread {resources.cpus} -I {input.bam} -R {input.ref} \
            -L 3R --out {output.vcf} --genotype_likelihoods_model BOTH \
            --genotyping_mode DISCOVERY --heterozygosity 0.015 --heterozygosity_stdev 0.05 --indel_heterozygosity 0.001 --downsampling_type BY_SAMPLE \
            -dcov 250 --output_mode EMIT_ALL_SITES --min_base_quality_score 17 -stand_call_conf 0.0 -contamination 0.0 -A DepthPerAlleleBySample \
            -A RMSMappingQuality -A Coverage -A FisherStrand -A StrandOddsRatio -A BaseQualityRankSumTest -A MappingQualityRankSumTest -A QualByDepth \
            -A ReadPosRankSumTest -XA ExcessHet -XA InbreedingCoeff -XA MappingQualityZero -XA HaplotypeScore -XA SpanningDeletions -XA ChromosomeCounts \
            &> {log}
        """

###############################################################################
rule bgzip_vcfs:
    input:
        rules.unifiedgenotyper.output.vcf
    output:
        gzip = "results/05_Variants/variants_X.vcf.gz",
    resources: cpus=1, mem_mb=4000, time_min=1200
    params: partition = 'fast'
    log:
        "results/11_Reports/bgzip/variants_X.vcf.gz.log",
    shell:
        config["MODULES"]["HTSLIB"]+"""
            bgzip -c --threads {threads} {input} > {output.gzip} {log}
        """

###############################################################################
rule bgzip_vcfs_2L:
    input:
        rules.unifiedgenotyper_2L.output.vcf
    output:
        gzip = "results/05_Variants/variants_2L.vcf.gz",
    resources: cpus=1, mem_mb=4000, time_min=1200
    params: partition = 'fast'
    log:
        "results/11_Reports/bgzip/variants_2L.vcf.gz.log",
    shell:
        config["MODULES"]["HTSLIB"]+"""
            bgzip -c --threads {threads} {input} > {output.gzip} {log}
        """

###############################################################################
rule bgzip_vcfs_2R:
    input:
        rules.unifiedgenotyper_2R.output.vcf
    output:
        gzip = "results/05_Variants/variants_2R.vcf.gz",
    resources: cpus=1, mem_mb=4000, time_min=1200
    params: partition = 'fast'
    log:
        "results/11_Reports/bgzip/variants_2R.vcf.gz.log",
    shell:
        config["MODULES"]["HTSLIB"]+"""
            bgzip -c --threads {threads} {input} > {output.gzip} {log}
        """

###############################################################################
rule bgzip_vcfs_3L:
    input:
        rules.unifiedgenotyper_3L.output.vcf
    output:
        gzip = "results/05_Variants/variants_3L.vcf.gz",
    resources: cpus=1, mem_mb=4000, time_min=1200
    params: partition = 'fast'
    log:
        "results/11_Reports/bgzip/variants_3L.vcf.gz.log",
    shell:
        config["MODULES"]["HTSLIB"]+"""
            bgzip -c --threads {threads} {input} > {output.gzip} {log}
        """

###############################################################################
rule bgzip_vcfs_3R:
    input:
        rules.unifiedgenotyper_3R.output.vcf
    output:
        gzip = "results/05_Variants/variants_3R.vcf.gz",
    resources: cpus=1, mem_mb=4000, time_min=1200
    params: partition = 'fast'
    log:
        "results/11_Reports/bgzip/variants_3R.vcf.gz.log",
    shell:
        config["MODULES"]["HTSLIB"]+"""
            bgzip -c --threads {threads} {input} > {output.gzip} {log}
        """

###############################################################################
rule indexfeaturefile:
    #threads: 4    
    message:
        "Indexing VCFs file for VariantToTable"
    input:
        vcf = rules.bgzip_vcfs.output.gzip
    output:
        indexvcf = "results/05_Variants/variants_X.vcf.gz.tbi",
    resources: cpus=1, mem_mb=4000, time_min=120
    params: partition = 'fast'
    log:
        "results/11_Reports/indexfeaturefile/variants_X.vcf.gz.tbi.log",
    shell:
        config["MODULES"]["GATK4"]+"""
            gatk IndexFeatureFile -I {input.vcf} -O {output.indexvcf} > {log} 2>&1 || true
        """

###############################################################################
rule indexfeaturefile_2L:
    #threads: 4    
    message:
        "Indexing VCFs file for VariantToTable"
    input:
        vcf = rules.bgzip_vcfs_2L.output.gzip
    output:
        indexvcf = "results/05_Variants/variants_2L.vcf.gz.tbi",
    resources: cpus=1, mem_mb=4000, time_min=120
    params: partition = 'fast'
    log:
        "results/11_Reports/indexfeaturefile/variants_2L.vcf.gz.tbi.log",
    shell:
        config["MODULES"]["GATK4"]+"""
            gatk IndexFeatureFile -I {input.vcf} -O {output.indexvcf} > {log} 2>&1 || true
        """

###############################################################################
rule indexfeaturefile_2R:
    #threads: 4    
    message:
        "Indexing VCFs file for VariantToTable"
    input:
        vcf = rules.bgzip_vcfs_2R.output.gzip
    output:
        indexvcf = "results/05_Variants/variants_2R.vcf.gz.tbi",
    resources: cpus=1, mem_mb=4000, time_min=120
    params: partition = 'fast'
    log:
        "results/11_Reports/indexfeaturefile/variants_2R.vcf.gz.tbi.log",
    shell:
        config["MODULES"]["GATK4"]+"""
            gatk IndexFeatureFile -I {input.vcf} -O {output.indexvcf} > {log} 2>&1 || true
        """

###############################################################################
rule indexfeaturefile_3L:
    #threads: 4    
    message:
        "Indexing VCFs file for VariantToTable"
    input:
        vcf = rules.bgzip_vcfs_3L.output.gzip
    output:
        indexvcf = "results/05_Variants/variants_3L.vcf.gz.tbi",
    resources: cpus=1, mem_mb=4000, time_min=120
    params: partition = 'fast'
    log:
        "results/11_Reports/indexfeaturefile/variants_3L.vcf.gz.tbi.log",
    shell:
        config["MODULES"]["GATK4"]+"""
            gatk IndexFeatureFile -I {input.vcf} -O {output.indexvcf} > {log} 2>&1 || true
        """

###############################################################################
rule indexfeaturefile_3R:
    #threads: 4    
    message:
        "Indexing VCFs file for VariantToTable"
    input:
        vcf = rules.bgzip_vcfs_3R.output.gzip
    output:
        indexvcf = "results/05_Variants/variants_3R.vcf.gz.tbi",
    resources: cpus=1, mem_mb=4000, time_min=120
    params: partition = 'fast'
    log:
        "results/11_Reports/indexfeaturefile/variants_3R.vcf.gz.tbi.log",
    shell:
        config["MODULES"]["GATK4"]+"""
            gatk IndexFeatureFile -I {input.vcf} -O {output.indexvcf} > {log} 2>&1 || true
        """

###############################################################################
def get_vcf_list(list):
    vcf_list = " --variant ".join(list)
    return(f"--variant {vcf_list}")


###############################################################################
rule variantstotable:
    # Aim:  Call variants in sequence data. The following parameters comes from the MalariaGEN
    # Use:  gatk VariantsToTable\
    #       -V input.vcf \
    #       -F CHROM -F POS -F TYPE -GF AD \
    #       -O output.table
    #threads: 2
    message:
        "VariantsToTable for chromosome X"
    input:
        vcf = rules.bgzip_vcfs.output.gzip,
        indexvcf = rules.indexfeaturefile.output.indexvcf,
    output:
        table = "results/05_Variants/vcf_X.table",
    resources: cpus=1, mem_mb=4000, time_min=120
    params: partition = 'long'
    log:
        "results/11_Reports/variantstotable/vcf_X.log",
    shell:
        config["MODULES"]["GATK4"]+"""
            gatk VariantsToTable -V {input.vcf} -O {output.table} -F CHROM -F POS -F QUAL -F DP -F FS -F MQ -F MQRankSum -F ReadPosRankSum -F QD -F SOR &> {log}
        """

###############################################################################
rule variantstotable_2L:
    # Aim:  Call variants in sequence data. The following parameters comes from the MalariaGEN
    # Use:  gatk VariantsToTable\
    #       -V input.vcf \
    #       -F CHROM -F POS -F TYPE -GF AD \
    #       -O output.table
    #threads: 2
    message:
        "VariantsToTable for chromosome 2L"
    input:
        vcf = rules.bgzip_vcfs_2L.output.gzip,
        indexvcf = rules.indexfeaturefile_2L.output.indexvcf,
    output:
        table = "results/05_Variants/vcf_2L.table",
    resources: cpus=1, mem_mb=4000, time_min=120
    params: partition = 'long'
    log:
        "results/11_Reports/variantstotable/vcf_2L.log",
    shell:
        config["MODULES"]["GATK4"]+"""
            gatk VariantsToTable -V {input.vcf} -O {output.table} -F CHROM -F POS -F QUAL -F DP -F FS -F MQ -F MQRankSum -F ReadPosRankSum -F QD -F SOR &> {log}
        """

###############################################################################
rule variantstotable_2R:
    # Aim:  Call variants in sequence data. The following parameters comes from the MalariaGEN
    # Use:  gatk VariantsToTable\
    #       -V input.vcf \
    #       -F CHROM -F POS -F TYPE -GF AD \
    #       -O output.table
    #threads: 2
    message:
        "VariantsToTable for chromosome 2R"
    input:
        vcf = rules.bgzip_vcfs_2R.output.gzip,
        indexvcf = rules.indexfeaturefile_2R.output.indexvcf,
    output:
        table = "results/05_Variants/vcf_2R.table",
    resources: cpus=1, mem_mb=4000, time_min=120
    params: partition = 'long'
    log:
        "results/11_Reports/variantstotable/vcf_2R.log",
    shell:
        config["MODULES"]["GATK4"]+"""
            gatk VariantsToTable -V {input.vcf} -O {output.table} -F CHROM -F POS -F QUAL -F DP -F FS -F MQ -F MQRankSum -F ReadPosRankSum -F QD -F SOR &> {log}
        """

###############################################################################
rule variantstotable_3L:
    # Aim:  Call variants in sequence data. The following parameters comes from the MalariaGEN
    # Use:  gatk VariantsToTable\
    #       -V input.vcf \
    #       -F CHROM -F POS -F TYPE -GF AD \
    #       -O output.table
    #threads: 2
    message:
        "VariantsToTable for chromosome 3L"
    input:
        vcf = rules.bgzip_vcfs_3L.output.gzip,
        indexvcf = rules.indexfeaturefile_3L.output.indexvcf,
    output:
        table = "results/05_Variants/vcf_3L.table",
    resources: cpus=1, mem_mb=4000, time_min=120
    params: partition = 'long'
    log:
        "results/11_Reports/variantstotable/vcf_3L.log",
    shell:
        config["MODULES"]["GATK4"]+"""
            gatk VariantsToTable -V {input.vcf} -O {output.table} -F CHROM -F POS -F QUAL -F DP -F FS -F MQ -F MQRankSum -F ReadPosRankSum -F QD -F SOR &> {log}
        """

###############################################################################
rule variantstotable_3R:
    # Aim:  Call variants in sequence data. The following parameters comes from the MalariaGEN
    # Use:  gatk VariantsToTable\
    #       -V input.vcf \
    #       -F CHROM -F POS -F TYPE -GF AD \
    #       -O output.table
    #threads: 2
    message:
        "VariantsToTable for chromosome 3R"
    input:
        vcf = rules.bgzip_vcfs_3R.output.gzip,
        indexvcf = rules.indexfeaturefile_3R.output.indexvcf,
    output:
        table = "results/05_Variants/vcf_3R.table",
    resources: cpus=1, mem_mb=4000, time_min=120
    params: partition = 'long'
    log:
        "results/11_Reports/variantstotable/vcf_3R.log",
    shell:
        config["MODULES"]["GATK4"]+"""
            gatk VariantsToTable -V {input.vcf} -O {output.table} -F CHROM -F POS -F QUAL -F DP -F FS -F MQ -F MQRankSum -F ReadPosRankSum -F QD -F SOR &> {log}
        """

###############################################################################
rule vcf_stats_X:
    message: "Execute rule VCF stats for chromosome X"
    input:
        vcf_file_X = rules.bgzip_vcfs.output.gzip
    output:
        freq = 'results/06_SNP_calling_stats/vcf_X.frq',
        depth = 'results/06_SNP_calling_stats/vcf_X.idepth',
        depth_mean = 'results/06_SNP_calling_stats/vcf_X.ldepth.mean',
        qual = 'results/06_SNP_calling_stats/vcf_X.lqual',
        missing_ind = 'results/06_SNP_calling_stats/vcf_X.imiss',
        miss = 'results/06_SNP_calling_stats/vcf_X.lmiss',
    resources: cpus=1, mem_mb=4000, time_min=1200
    params: 
        partition = 'fast',
        outdir = 'results/'
    log:
            error =  'results/11_Reports/vcf_stats/vcftools_vcf_X.e',
            output = 'results/11_Reports/vcf_stats/vcftools_vcf_X.o'        
    shell:
        config["MODULES"]["VCFTOOLS"]+" && "+config["MODULES"]["SAMTOOLS"]+"""
            vcftools --gzvcf {input.vcf_file_X}  --remove-indels --freq2 --max-alleles 2 --stdout 1> {output.freq}
            vcftools --gzvcf {input.vcf_file_X}  --remove-indels --depth --stdout 1> {output.depth}
            vcftools --gzvcf {input.vcf_file_X}  --remove-indels --site-mean-depth --stdout 1> {output.depth_mean}
            vcftools --gzvcf {input.vcf_file_X}  --remove-indels --site-quality --stdout 1> {output.qual}
            vcftools --gzvcf {input.vcf_file_X}  --remove-indels --missing-indv --stdout 1> {output.missing_ind}
            vcftools --gzvcf {input.vcf_file_X}  --remove-indels --missing-site --stdout 1> {output.miss}
        """

###############################################################################
rule vcf_stats_2L:
    message: "Execute rule VCF stats for chromosome 2L"
    input:
        vcf_file_2L = rules.bgzip_vcfs_2L.output.gzip
    output:
        freq = 'results/06_SNP_calling_stats/vcf_2L.frq',
        depth = 'results/06_SNP_calling_stats/vcf_2L.idepth',
        depth_mean = 'results/06_SNP_calling_stats/vcf_2L.ldepth.mean',
        qual = 'results/06_SNP_calling_stats/vcf_2L.lqual',
        missing_ind = 'results/06_SNP_calling_stats/vcf_2L.imiss',
        miss = 'results/06_SNP_calling_stats/vcf_2L.lmiss',
    resources: cpus=1, mem_mb=4000, time_min=1200
    params: 
        partition = 'fast',
        outdir = 'results/'
    log:
            error =  'results/11_Reports/vcf_stats/vcftools_vcf_2L.e',
            output = 'results/11_Reports/vcf_stats/vcftools_vcf_2L.o'        
    shell:
        config["MODULES"]["VCFTOOLS"]+" && "+config["MODULES"]["SAMTOOLS"]+"""
            vcftools --gzvcf {input.vcf_file_2L}  --remove-indels --freq2 --max-alleles 2 --stdout 1> {output.freq}
            vcftools --gzvcf {input.vcf_file_2L}  --remove-indels --depth --stdout 1> {output.depth}
            vcftools --gzvcf {input.vcf_file_2L}  --remove-indels --site-mean-depth --stdout 1> {output.depth_mean}
            vcftools --gzvcf {input.vcf_file_2L}  --remove-indels --site-quality --stdout 1> {output.qual}
            vcftools --gzvcf {input.vcf_file_2L}  --remove-indels --missing-indv --stdout 1> {output.missing_ind}
            vcftools --gzvcf {input.vcf_file_2L}  --remove-indels --missing-site --stdout 1> {output.miss}
        """

###############################################################################
rule vcf_stats_2R:
    message: "Execute rule VCF stats for chromosome 2R"
    input:
        vcf_file_2R = rules.bgzip_vcfs_2R.output.gzip
    output:
        freq = 'results/06_SNP_calling_stats/vcf_2R.frq',
        depth = 'results/06_SNP_calling_stats/vcf_2R.idepth',
        depth_mean = 'results/06_SNP_calling_stats/vcf_2R.ldepth.mean',
        qual = 'results/06_SNP_calling_stats/vcf_2R.lqual',
        missing_ind = 'results/06_SNP_calling_stats/vcf_2R.imiss',
        miss = 'results/06_SNP_calling_stats/vcf_2R.lmiss',
    resources: cpus=1, mem_mb=4000, time_min=1200
    params:
        partition = 'fast',
        outdir = 'results/'
    log:
            error =  'results/11_Reports/vcf_stats/vcftools_vcf_2R.e',
            output = 'results/11_Reports/vcf_stats/vcftools_vcf_2R.o'        
    shell:
        config["MODULES"]["VCFTOOLS"]+" && "+config["MODULES"]["SAMTOOLS"]+"""
            vcftools --gzvcf {input.vcf_file_2R}  --remove-indels --freq2 --max-alleles 2 --stdout 1> {output.freq}
            vcftools --gzvcf {input.vcf_file_2R}  --remove-indels --depth --stdout 1> {output.depth}
            vcftools --gzvcf {input.vcf_file_2R}  --remove-indels --site-mean-depth --stdout 1> {output.depth_mean}
            vcftools --gzvcf {input.vcf_file_2R}  --remove-indels --site-quality --stdout 1> {output.qual}
            vcftools --gzvcf {input.vcf_file_2R}  --remove-indels --missing-indv --stdout 1> {output.missing_ind}
            vcftools --gzvcf {input.vcf_file_2R}  --remove-indels --missing-site --stdout 1> {output.miss}
        """

###############################################################################
rule vcf_stats_3L:
    message: "Execute rule VCF stats for chromosome 2L"
    input:
        vcf_file_3L = rules.bgzip_vcfs_3L.output.gzip
    output:
        freq = 'results/06_SNP_calling_stats/vcf_3L.frq',
        depth = 'results/06_SNP_calling_stats/vcf_3L.idepth',
        depth_mean = 'results/06_SNP_calling_stats/vcf_3L.ldepth.mean',
        qual = 'results/06_SNP_calling_stats/vcf_3L.lqual',
        missing_ind = 'results/06_SNP_calling_stats/vcf_3L.imiss',
        miss = 'results/06_SNP_calling_stats/vcf_3L.lmiss',
    resources: cpus=1, mem_mb=4000, time_min=1200
    params:
        partition = 'fast',
        outdir = 'results/'
    log:
            error =  'results/11_Reports/vcf_stats/vcftools_vcf_3L.e',
            output = 'results/11_Reports/vcf_stats/vcftools_vcf_3L.o'        
    shell:
        config["MODULES"]["VCFTOOLS"]+" && "+config["MODULES"]["SAMTOOLS"]+"""
            vcftools --gzvcf {input.vcf_file_3L}  --remove-indels --freq2 --max-alleles 2 --stdout 1> {output.freq}
            vcftools --gzvcf {input.vcf_file_3L}  --remove-indels --depth --stdout 1> {output.depth}
            vcftools --gzvcf {input.vcf_file_3L}  --remove-indels --site-mean-depth --stdout 1> {output.depth_mean}
            vcftools --gzvcf {input.vcf_file_3L}  --remove-indels --site-quality --stdout 1> {output.qual}
            vcftools --gzvcf {input.vcf_file_3L}  --remove-indels --missing-indv --stdout 1> {output.missing_ind}
            vcftools --gzvcf {input.vcf_file_3L}  --remove-indels --missing-site --stdout 1> {output.miss}
        """

###############################################################################
rule vcf_stats_3R:
    message: "Execute rule VCF stats for chromosome 3R"
    input:
        vcf_file_3R = rules.bgzip_vcfs_3R.output.gzip
    output:
        freq = 'results/06_SNP_calling_stats/vcf_3R.frq',
        depth = 'results/06_SNP_calling_stats/vcf_3R.idepth',
        depth_mean = 'results/06_SNP_calling_stats/vcf_3R.ldepth.mean',
        qual = 'results/06_SNP_calling_stats/vcf_3R.lqual',
        missing_ind = 'results/06_SNP_calling_stats/vcf_3R.imiss',
        miss = 'results/06_SNP_calling_stats/vcf_3R.lmiss',
    resources: cpus=1, mem_mb=4000, time_min=1200
    params:
        partition = 'fast',
        outdir = 'results/'
    log:
            error =  'results/11_Reports/vcf_stats/vcftools_vcf_3R.e',
            output = 'results/11_Reports/vcf_stats/vcftools_vcf_3R.o'        
    shell:
        config["MODULES"]["VCFTOOLS"]+" && "+config["MODULES"]["SAMTOOLS"]+"""
            vcftools --gzvcf {input.vcf_file_3R}  --remove-indels --freq2 --max-alleles 2 --stdout 1> {output.freq}
            vcftools --gzvcf {input.vcf_file_3R}  --remove-indels --depth --stdout 1> {output.depth}
            vcftools --gzvcf {input.vcf_file_3R}  --remove-indels --site-mean-depth --stdout 1> {output.depth_mean}
            vcftools --gzvcf {input.vcf_file_3R}  --remove-indels --site-quality --stdout 1> {output.qual}
            vcftools --gzvcf {input.vcf_file_3R}  --remove-indels --missing-indv --stdout 1> {output.missing_ind}
            vcftools --gzvcf {input.vcf_file_3R}  --remove-indels --missing-site --stdout 1> {output.miss}
        """

###############################################################################
rule report_vcf:
    message: "Execute rule report_vcf for chromosome X",
    input:
        freq = rules.vcf_stats_X.output.freq,
        depth = rules.vcf_stats_X.output.depth,
        depth_mean = rules.vcf_stats_X.output.depth_mean,
        qual = rules.vcf_stats_X.output.qual,
        missing_ind = rules.vcf_stats_X.output.missing_ind,
        miss = rules.vcf_stats_X.output.miss
    params:
        outdir = 'results/',
        partition = 'fast'
    output:
        report = "results/report_vcf_X.html",
    resources: cpus=1, mem_mb=24000, time_min=120
    log:
        error =  'results/11_Reports/report_vcf/report_X.e',
        output = 'results/11_Reports/report_vcf/report_X.o',
    script:
        """scripts/report_vcf.Rmd"""

###############################################################################
rule report_vcf_2L:
    message: "Execute rule report_vcf for chromosome 2L",
    input:
        freq = rules.vcf_stats_2L.output.freq,
        depth = rules.vcf_stats_2L.output.depth,
        depth_mean = rules.vcf_stats_2L.output.depth_mean,
        qual = rules.vcf_stats_2L.output.qual,
        missing_ind = rules.vcf_stats_2L.output.missing_ind,
        miss = rules.vcf_stats_2L.output.miss
    params:
        outdir = 'results/',
        partition = 'fast'
    output:
        report = "results/report_vcf_2L.html",
    resources: cpus=1, mem_mb=24000, time_min=120
    log:
        error =  'results/11_Reports/report_vcf/report_2L.e',
        output = 'results/11_Reports/report_vcf/report_2L.o'
    script:
        """scripts/report_vcf.Rmd"""

###############################################################################
rule report_vcf_2R:
    message: "Execute rule report_vcf for chromosome 2R",
    input:
        freq = rules.vcf_stats_2R.output.freq,
        depth = rules.vcf_stats_2R.output.depth,
        depth_mean = rules.vcf_stats_2R.output.depth_mean,
        qual = rules.vcf_stats_2R.output.qual,
        missing_ind = rules.vcf_stats_2R.output.missing_ind,
        miss = rules.vcf_stats_2R.output.miss
    params:
        outdir = 'results/',
        partition = 'fast'
    output:
        report = "results/report_vcf_2R.html",
    resources: cpus=1, mem_mb=24000, time_min=120
    log:
        error =  'results/11_Reports/report_vcf/report_2R.e',
        output = 'results/11_Reports/report_vcf/report_2R.o'
    script:
        """scripts/report_vcf.Rmd"""

###############################################################################
rule report_vcf_3L:
    message: "Execute rule report_vcf for chromosome 3L",
    input:
        freq = rules.vcf_stats_3L.output.freq,
        depth = rules.vcf_stats_3L.output.depth,
        depth_mean = rules.vcf_stats_3L.output.depth_mean,
        qual = rules.vcf_stats_3L.output.qual,
        missing_ind = rules.vcf_stats_3L.output.missing_ind,
        miss = rules.vcf_stats_3L.output.miss
    params:
        outdir = 'results/',
        partition = 'fast'
    output:
        report = "results/report_vcf_3L.html",
    resources: cpus=1, mem_mb=24000, time_min=120
    log:
        error =  'results/11_Reports/report_vcf/report_3L.e',
        output = 'results/11_Reports/report_vcf/report_3L.o'
    script:
        """scripts/report_vcf.Rmd"""

###############################################################################
rule report_vcf_3R:
    message: "Execute rule report_vcf for chromosome 3R",
    input:
        freq = rules.vcf_stats_3R.output.freq,
        depth = rules.vcf_stats_3R.output.depth,
        depth_mean = rules.vcf_stats_3R.output.depth_mean,
        qual = rules.vcf_stats_3R.output.qual,
        missing_ind = rules.vcf_stats_3R.output.missing_ind,
        miss = rules.vcf_stats_3R.output.miss
    params:
        outdir = 'results/',
        partition = 'fast'
    output:
        report = "results/report_vcf_3R.html",
    resources: cpus=1, mem_mb=24000, time_min=120
    log:
        error =  'results/11_Reports/report_vcf/report_3R.e',
        output = 'results/11_Reports/report_vcf/report_3R.o'
    script:
        """scripts/report_vcf.Rmd"""



