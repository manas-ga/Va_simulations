

################################################################################################################
### Loop through output files, extracting mutations and genomes using the extract_mutations_genomes function ###
################################################################################################################

# Loop through all the generations where SLiM generated outputs, generate appropriate filenames for SLiM output files
# For each file, extract mutations and genomes
# Enter the starting generation, the ending generation and the frequency of the output

start_gen = 101000
end_gen = 200000
output_freq = 1000


# Import the function to extract mutations and genomes from SLiM outputs (currently stored in the "a_Functions" directory)

from a_Functions.Function_to_extract_mutations_and_genomes_from_SLiM_output import extract_mutations_genomes


for gen in range(start_gen, end_gen + 1, output_freq):
    file = f"Output_{gen}.txt"
    extract_mutations_genomes(slim_output_path = f"b_Interim_files/SLiM_Outputs/{file}", mutations_path = f"b_Interim_files/Mutations/mutations_{gen}.txt", all_genomes_path = f"b_Interim_files/Genomes/all_genomes_{gen}.txt", c_matrix_path = f"b_Interim_files/C_Matrices/c_matrix_{gen}.csv")
    print(f"{int((gen-start_gen)/(end_gen-start_gen)*100)} % complete")
    




