#!/bin/bash

#SBATCH --job-name=lost_vA
#SBATCH --output=/mnt/hel/obbard/Va_simulations/analyses/b_Interim_files/std_out/job_%A_%a.log   # Store logs in a custom directory
#SBATCH --open-mode=append                                                                       # Append output if the file already exists
#SBATCH --array=1-80%10                                                                       # Run replicate tasks
#SBATCH --ntasks=10
#SBATCH --mem=30G

# Path to Conda (if Conda isn't initialized by default)
CONDA_PATH="/home/msamant/miniconda3"  # Update to where Miniconda/Anaconda is installed
source $CONDA_PATH/etc/profile.d/conda.sh  # Initialize Conda in the job script

# Activate the Conda environment
conda activate marun

Param=`cat lost_va_IDs.txt | awk "NR==${SLURM_ARRAY_TASK_ID}"`

Rscript lost_va_script.R $Param
