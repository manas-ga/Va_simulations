#!/bin/bash

#SBATCH --job-name=Set_9
#SBATCH --output=/mnt/hel/obbard/Va_simulations/analyses/b_Interim_files/std_out/job_%A_%a.log   # Store logs in a custom directory
#SBATCH --open-mode=append                                                                       # Append output if the file already exists
#SBATCH --array=1-100%5                                                                          # Run replicate tasks
#SBATCH --ntasks=12
#SBATCH --mem=250G
#SBATCH --nodelist=ac3-n6                                                                        # Only submit to n6 (has large RAM)

# Path to Conda (if Conda isn't initialized by default)
CONDA_PATH="/home/msamant/miniconda3"  # Update to where Miniconda/Anaconda is installed
source $CONDA_PATH/etc/profile.d/conda.sh  # Initialize Conda in the job script

# Activate the Conda environment
conda activate marun

Param=`cat analysis_param_grid.txt | awk "NR==${SLURM_ARRAY_TASK_ID}"`

Rscript Analysis_script.R $Param
