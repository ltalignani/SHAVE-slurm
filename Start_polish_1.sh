#!/bin/bash

###################configuration slurm##############################
# SHAVE
#sbatch --account=aedes_amplicon
#SBATCH --job-name=SHAVE_POLISHING
#SBATCH --output=shave_polishing.%A.out
#SBATCH --error=shave_polishing.%A.err
#SBATCH --partition=long
#SBATCH --nodes=1
#sbatch --cpus-per-task=4
#SBATCH --mem=4GB
# Define email for script execution
#SBATCH --mail-user=loic.talignani@ird.fr
# Define type notifications (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-type=ALL
###################################################################

# USAGE: sbatch Start_polishing.sh

###Charge module
module purge
module load samtools/1.15.1 snakemake/7.7.0 slurm-drmaa/1.0.8 picard/2.23.5 graphviz/2.40.1
echo $DRMAA_LIBRARY_PATH

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
echo -e "${blue}Name${nc} __________________ Start_polishing.sh"
echo -e "${blue}Author${nc} ________________ Loïc Talignani"
echo -e "${blue}Affiliation${nc} ___________ UMR_MIVEGEC"
echo -e "${blue}Aim${nc} ___________________ Bash script for ${red}SH${nc}ort-read ${red}A${nc}lignment pipeline for ${red}VE${nc}ctor v.1"
echo -e "${blue}Date${nc} __________________ 2023.01.25"
echo -e "${blue}Run${nc} ___________________ bash Start_polishing.sh"
echo -e "${blue}Latest Modification${nc} ___ Updated for IFB-core cluster"

# Operating System
os="Linux"

# Set working directory
workdir=$(pwd)            #$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
max_threads="24"

echo -e "Workdir is "${workdir}

###### Call snakemake pipeline ######
echo ""
echo -e "${green}------------------------------------------------------------------------${nc}"
echo -e "${green}###########${nc} ${red}SNAKEMAKE PIPELINE START${nc} ${green}###########${nc}"
echo -e "${green}------------------------------------------------------------------------${nc}"
echo ""

# slurm directive by rule
CLUSTER_CONFIG=cluster.yaml
# sbatch directive to pass to snakemake
CLUSTER='sbatch --account=aedes_amplicon --partition={cluster.partition} --mem={cluster.mem} -c {cluster.cpus} -o {cluster.output} -e {cluster.error}'
DRMAA='--account=aedes_amplicon --partition={cluster.partition} --mem={cluster.mem} -c {cluster.cpus} -o {cluster.output} -e {cluster.error}'
# Maximum number of jobs to be submitted at a time (see cluster limitation)
MAX_JOBS=250

# Clean up .snakemake directory
#rm -rf .snakemake

echo -e "${blue}Unlocking working directory:${nc}"
echo ""
# Specify working directory (relative paths in the snakefile will use this as their origin).
# The slurm profile for a cluster under a SLURM workload manager.
# The workflow definition in form of a snakefile.
# Set or overwrite values in the workflow config object.
# Re-run all jobs the output of which is recognized as incomplete.
# Remove a lock on the working directory.
snakemake \
    --directory ${workdir}/ \
    --snakefile workflow/rules/shave_polish_1.smk \
    --config os=${os} \
    --rerun-incomplete \
    --unlock \
    --use-conda \
    --conda-frontend conda # Default "mamba", recommended because much faster, but : "Library not loaded: @rpath/libarchive.13.dylib"

echo ""
echo -e "${blue}Conda environments list:${nc}"
echo ""
# Specify working directory (relative paths in the snakefile will use this as their origin).
# The slurm profile for a cluster under a SLURM workload manager.
# The workflow definition in form of a snakefile.
# Set or overwrite values in the workflow config object.
# Re-run all jobs the output of which is recognized as incomplete.
# List all conda environments and their location on disk.
snakemake \
    --directory ${workdir}/ \
    --snakefile workflow/rules/shave_polish_1.smk \
    --config os=${os} \
    --rerun-incomplete \
    --use-conda \
    --list-conda-envs \
    --conda-frontend conda

echo ""
echo -e "${blue}Conda environments update:${nc}"
echo ""
# Specify working directory (relative paths in the snakefile will use this as their origin).
# The slurm profile for a cluster under a SLURM workload manager.
# The workflow definition in form of a snakefile.
# Re-run all jobs the output of which is recognized as incomplete.
# Set or overwrite values in the workflow config object.
# Cleanup unused conda environments.
snakemake \
    --directory ${workdir}/ \
    --snakefile workflow/rules/shave_polish_1.smk \
    --cores ${max_threads} \
    --config os=${os} \
    --rerun-incomplete \
    --conda-cleanup-envs\
    --conda-frontend conda \
    --use-envmodules \
    --use-conda

# snakemake --directory ./ --snakefile shave.smk --config os="Linux" --rerun-incomplete --conda-cleanup-envs --use-conda --list-conda-envs --conda-frontend conda --use-envmodules
echo ""
echo -e "${blue}Conda environments setup:${nc}"
echo ""
# Specify working directory (relative paths in the snakefile will use this as their origin).
# The slurm profile for a cluster under a SLURM workload manager.
# The workflow definition in form of a snakefile.
# Set or overwrite values in the workflow config object.
# Re-run all jobs the output of which is recognized as incomplete.
# If defined in the rule, run job in a conda environment.
# If specified, only creates the job-specific conda environments then exits. The –use-conda flag must also be set.
# If mamba package manager is not available, or if you still prefer to use conda, you can enforce that with this setting (default: 'mamba').
snakemake \
    --directory ${workdir}/ \
    --snakefile workflow/rules/shave_polish_1.smk \
    --cores ${max_threads} \
    --config os=${os} \
    --rerun-incomplete \
    --use-envmodules \
    --conda-frontend conda
    --conda-create-envs-only \
    --use-conda

echo ""
echo -e "${blue}Let's run!${nc}"
echo ""
# Specify working directory (relative paths in the snakefile will use this as their origin).
# The slurm profile for a cluster under a SLURM workload manager.
# The workflow definition in form of a snakefile.
# Set or overwrite values in the workflow config object.
# Re-run all jobs the output of which is recognized as incomplete.
# Go on with independent jobs if a job fails.
# If defined in the rule, run job in a conda environment.
# Tell the scheduler to assign creation of given targets (and all their dependencies) highest priority.
# Print out the shell commands that will be executed.
snakemake \
    --directory ./ \
    --snakefile workflow/rules/shave_polish_1.smk \
    --jobs=250 \
    --cores ${max_threads} \
    --config os="Linux" \
    --rerun-incomplete \
    --keep-going \
    --conda-frontend conda \
    --printshellcmds \
    --drmaa " --account=aedes_amplicon --partition={cluster.partition} --mem={cluster.mem} -c {cluster.cpus} -o {cluster.output} -e {cluster.error}" \
    --cluster-config cluster.yaml \
    --configfile config.yaml
    --use-envmodules \
    --use-conda

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
	snakemake \
	    --directory ${workdir}/ \
        --cluster-config cluster.yaml \
        --snakefile workflow/rules/shave_polish_1.smk \
        --${graph} | \
	    dot -T${extention} > \
		${workdir}/results/10_Graphs/${graph}.${extention} ;
    done ;
done

snakemake \
    --directory ${workdir} \
    --cluster-config cluster.yaml \
    --snakefile workflow/rules/shave_polish_1.smk \
    --summary > ${workdir}/results/11_Reports/files_summary.txt

cp ${workdir}/config/config.yaml ${workdir}/results/11_Reports/config.yaml

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
