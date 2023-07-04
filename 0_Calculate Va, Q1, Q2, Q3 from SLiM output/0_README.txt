##########################################################################################################
################### Compute Va, Q1, Q2, and Q3 from the output of a W-F SLiM simulations #################
##########################################################################################################

## This folder contains three scripts:

### "1_W-F_sim.slim" ###

# This is a SLiM file that runs a Wright-Fisher simulation.
# At the end of the simulation a sample of 100 individuals is drawn from the simulated population
# The final output is stored in the "Interim_files" folder as "Output_WF_SLiM.txt"


### "2_Read mutations and genomes" ###

#This is a python file that reads "Output_WF_SLiM.txt" and generates three files that are stored in "Interim_files":
# 1. A text file of all mutations and their information in the entire population ("mutations.txt")
# 2. A text file of all genomes containing lists of their mutations in the sample drawn at the end of the simulation ("all_genomes.txt"). Note that this file is not used in downstream analyses.
# 3. A csv file containing the c matrix for genomes in the said sample ("c_matrix.csv"). Note that the genomes in the SLiM output are such that genomes 2*i and 2*i - 1 belong to individual i

### "3_Compute V1, Q1, Q2, Q3" ###

# This is an R file that reads "mutations.txt" and "c_matrix.csv", fits a multiple regression of allele counts on fitness, and calculates Va, Q1, Q2, and Q3