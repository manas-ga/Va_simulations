

################################################################################################################
### Loop through output files, extracting mutations and genomes using the extract_mutations_genomes function ###
################################################################################################################

# Loop through all the generations where SLiM generated outputs, generate appropriate filenames for SLiM output files
# For each file, extract mutations and genomes (from the sampled individuals, sample_size is parameter of the function) and store files in appropriate locations


# Enter the starting generation, the ending generation and the frequency of the output

start_gen = 1
end_gen = 5
output_freq = 1

# Enter workig directory

wd = "C:/Academics/Post-doc/Va_simulations/4_Va_theory/"

# Import the function to extract mutations and genomes from SLiM outputs and the function to extract just the mutations (currently stored in the "a_Functions" directory)

from a_Functions.Function_to_extract_mutations_and_genomes_from_SLiM_output import extract_mutations_genomes
from a_Functions.Function_to_extract_mutations import extract_mutations


# Extract just the mutations for all generations


for gen in range(start_gen, end_gen + 1, output_freq):
    file = f"Output_{gen}.txt"
    extract_mutations(slim_output_path = f"{wd}b_Interim_files/SLiM_Outputs/{file}", mutations_path = f"{wd}b_Interim_files/Mutations/mutations_{gen}.txt")
                      
    #print(f"{round((gen-(start_gen-output_freq))/(end_gen-(start_gen-output_freq))*100, 2)} % complete!")








#### The code below can be used if one wants the c matrix from the initial generation as well
# Otherwise ignore, as computing the c matrix takes time

# Extract mutations as well as the genomes for the first generation

#gen = start_gen

#file = f"Output_{gen}.txt"





# Extract mutations and genomes for all generations


for gen in range(start_gen, end_gen + 1, output_freq):
    file = f"Output_{gen}.txt"
    extract_mutations_genomes(slim_output_path = f"{wd}b_Interim_files/SLiM_Outputs/{file}",
                          mutations_path = f"{wd}b_Interim_files/Mutations/mutations_{gen}.txt",
                          all_genomes_path = f"{wd}b_Interim_files/Genomes/all_genomes_{gen}.txt",
                          c_matrix_path = f"{wd}b_Interim_files/C_Matrices/c_matrix_{gen}.csv",
                          sample_size = 1000)
                      
    print(f"{round((gen-(start_gen-output_freq))/(end_gen-(start_gen-output_freq))*100, 2)} % complete!")
    




