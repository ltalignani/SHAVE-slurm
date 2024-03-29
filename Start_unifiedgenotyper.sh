#!/bin/bash

###################configuration slurm##############################
# SHAVE
#sbatch --account=aedes_amplicon
#SBATCH --job-name=SHAVE
#SBATCH --output=shave.%A.out
#SBATCH --error=shave.%A.err
#SBATCH --partition=long
#SBATCH --nodes=1
#sbatch --cpus-per-task=1
#SBATCH --mem=1GB
# Define email for script execution
#SBATCH --mail-user=loic.talignani@ird.fr
# Define type notifications (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-type=ALL
###################################################################

# USAGE: sbatch Start_unifiedgenotyper.sh

###Charge module
module purge
module load snakemake/7.7.0
module load graphviz/2.40.1
module load r/4.2.1
#module load slurm-drmaa/1.0.8 
#echo $DRMAA_LIBRARY_PATH

# Activate private conda environment
#CONDA_ROOT=/shared/software/miniconda/
#source $CONDA_ROOT/etc/profile.d/conda.sh
#source ${HOME}/miniconda3/etc/profile.d/conda.sh

#conda activate gatk3

echo -e "${green} Current directory is $(pwd) ${nc}"

# set umask to avoid locking each other out of directories
umask 002

##### Colors ######
red="\033[1;31m"   # red
green="\033[1;32m" # green
ylo="\033[1;33m"   # yellow
blue="\033[1;34m"  # blue
nc="\033[0m"       # no color

###### About ######
echo ""
echo -e "${green}------------------------------------------------------------------------${nc}"
echo -e "${green}#####${nc} ${red}ABOUT${nc} ${green}#####${nc}"
echo -e "${green}-----------------${nc}"
echo ""
echo -e "${blue}Name${nc} __________________ Start_unifiedgenotyper.sh"
echo -e "${blue}Author${nc} ________________ Loïc Talignani"
echo -e "${blue}Affiliation${nc} ___________ UMR_MIVEGEC"
echo -e "${blue}Aim${nc} ___________________ Bash script for ${red}SH${nc}ort-read ${red}A${nc}lignment pipeline for ${red}VE${nc}ctor v.1"
echo -e "${blue}Date${nc} __________________ 2023.01.25"
echo -e "${blue}Run${nc} ___________________ bash Start_unifiedgenotyper.sh"
echo -e "${blue}Latest Modification${nc} ___ Updated for IFB-core cluster"

# Set working directory
workdir=$(pwd)            #$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
max_threads="30"

echo -e "Workdir is "${workdir}

###### Call snakemake pipeline ######
echo ""
echo -e "${green}------------------------------------------------------------------------${nc}"
echo -e "${green}###########${nc} ${red}SNAKEMAKE PIPELINE START${nc} ${green}###########${nc}"
echo -e "${green}------------------------------------------------------------------------${nc}"
echo ""

# slurm directive by rule
# CLUSTER_CONFIG=cluster.yaml
# # sbatch directive to pass to snakemake
# CLUSTER='sbatch --account=aedes_amplicon --partition={cluster.partition} --mem={cluster.mem} -c {cluster.cpus} -o {cluster.output} -e {cluster.error}'
# DRMAA='--account=aedes_amplicon --partition={cluster.partition} --mem={cluster.mem} -c {cluster.cpus} -o {cluster.output} -e {cluster.error}'
# # Maximum number of jobs to be submitted at a time (see cluster limitation)
# MAX_JOBS=250

# Clean up .snakemake directory
#rm -rf .snakemake

echo -e "${blue}Unlocking working directory:${nc}"
echo ""

snakemake --profile slurm/ --directory ${workdir}/ --snakefile workflow/rules/shave_unifiedgenotyper.smk --unlock

echo ""
echo -e "${blue}Let's run!${nc}"
echo ""

snakemake --profile slurm/ --snakefile workflow/rules/shave_unifiedgenotyper.smk --cores 30 --configfile config/config.yaml

# snakemake --profile slurm/ --snakefile workflow/rules/shave_unifiedgenotyper.smk --jobs=25 --cores 24 --rerun-incomplete --keep-going --conda-frontend conda --printshellcmds --configfile config.yaml --use-conda
#    --drmaa "--account=aedes_amplicon --partition={cluster.partition} --mem={cluster.mem} -c {cluster.cpus} -o {cluster.output} -e {cluster.error}" \
#    --cluster-config cluster.yaml \
# snakemake --directory ./ --snakefile shave.smk --jobs=250 --cores 24 --config os="Linux" --rerun-incomplete ---keep-going -use-conda --conda-frontend conda -prioritize multiqc --printshellcmds --cluster-config cluster.yaml --configfile config.yaml --use-envmodules --drmaa " --account=aedes_amplicon --partition={cluster.partition} --mem={cluster.mem} -c {cluster.cpus} -o {cluster.output} -e {cluster.error}"
# --cluster "sbatch --account=aedes_amplicon --partition={cluster.partition} --mem={cluster.mem} -c {cluster.cpus} -o {cluster.output} -e {cluster.error}" \

###### Create usefull graphs, summary and logs ######
echo ""
echo -e "${green}------------------------------------------------------------------------${nc}"
echo -e "${green}###########${nc} ${red}SNAKEMAKE PIPELINE LOGS${nc} ${green}############${nc}"
echo -e "${green}------------------------------------------------------------------------${nc}"
echo ""

mkdir ${workdir}/results/10_Graphs/ 2> /dev/null

graph_list="dag rulegraph filegraph"
extention_list="pdf png"

for graph in ${graph_list} ; do
    for extention in ${extention_list} ; do
	snakemake --profile slurm/ --snakefile workflow/rules/shave_unifiedgenotyper.smk --${graph} | dot -T${extention} > ${workdir}/results/10_Graphs/${graph}.${extention} ;
    done ;
done

snakemake --directory ${workdir} --profile slurm/ --snakefile workflow/rules/shave_unifiedgenotyper.smk --summary > ${workdir}/results/11_Reports/files_summary.txt

###### End managment ######
echo ""
echo -e "${blue}------------------------------------------------------------------------${nc}"
echo -e "${blue}##################${nc} ${red}SCRIPT END${nc} ${blue}###################${nc}"
echo -e "${blue}------------------------------------------------------------------------${nc}"
echo ""

find ${workdir}/results/ -type f -empty -delete                 # Remove empty file (like empty log)
find ${workdir}/results/ -type d -empty -delete                 # Remove empty directory

time_stamp_end=$(date +"%Y-%m-%d %H:%M")                        # Get date / hour ending analyzes
elapsed_time=${SECONDS}                                         # Get SECONDS counter
minutes=$((${elapsed_time}/60))                                 # / 60 = minutes
seconds=$((${elapsed_time}%60))                                    # % 60 = seconds

echo -e "${red}End Time${nc} ______________ ${time_stamp_end}"                                                             # Print analyzes ending time
echo -e "${red}Processing Time${nc} _______ ${ylo}${minutes}${nc} minutes and ${ylo}${seconds}${nc} seconds elapsed"       # Print total time elapsed

echo ""
echo -e "${blue}------------------------------------------------------------------------${nc}"
echo ""
