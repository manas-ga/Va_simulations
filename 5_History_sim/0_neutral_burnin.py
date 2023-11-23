## Modified from: https://tskit.dev/pyslim/docs/stable/vignette_coalescent_diversity.html


######## Python script that uses msprime and pyslim to generate a neutral burnin simulation ########

import pyslim
import msprime
import numpy
import sys

############## Output path #################

output_path = sys.argv[8] ### The output path is a command line argument to be fed by the controlling R script

###### Simulate a neutral ancestry using msprime ######

demog_model = msprime.Demography()
demog_model.add_population(initial_size=float(sys.argv[1]))

ots = msprime.sim_ancestry(samples=float(sys.argv[2]), sequence_length = float(sys.argv[3]), demography=demog_model, recombination_rate=float(sys.argv[4]))

# Add annotations for SLiM

ots = pyslim.annotate(ots, model_type="WF", tick=1, stage="late")

####### Add mutations #######

# Specifies that these mutations are going to be drawn using the infinite alleles model and are going to be assigned type "m2" in SLiM

mut_model = msprime.SLiMMutationModel(type=2)

ots = msprime.sim_mutations(ots, rate = float(sys.argv[5]), model = mut_model, keep = True)

print(f"The tree sequence now has {ots.num_mutations} mutations, at "
      f"{ots.num_sites} distinct sites.")


####### Assign selection coefficients #########

tables = ots.tables
tables.mutations.clear()

# Decide the relative frequencies of neutral vs selected mutations


freq_deleterious = 1000
freq_beneficial = int(freq_deleterious*float(sys.argv[9]))
freq_neutral = 0

# Creae a list of 0s and 1s in frequencies of neutral and selected mutations
# This list will be used to generate a random 0 or 1 to decide whether the mutation becomes neutral or selected

freq_array = [0]*freq_neutral + [-1]*freq_deleterious + [1]*freq_beneficial 

mut_map = {}
for m in ots.mutations():
  md_list = m.metadata["mutation_list"]
  slim_ids = m.derived_state.split(",")
  assert len(slim_ids) == len(md_list)
  for sid, md in zip(slim_ids, md_list):
    
    if sid not in mut_map:
      
      if sys.argv[10] == "g":
        
        # DFE is gamma
        mut_map[sid] = freq_array[numpy.random.randint(0, len(freq_array))]*(min(numpy.random.gamma(shape = float(sys.argv[6]), scale = float(sys.argv[7])), 1.0)) # The max() function prevents negative fitnesses

      if sys.argv[10] == "n":
        mut_map[sid] = numpy.random.normal(float(sys.argv[11]), float(sys.argv[12]))
        

      # DFE is just + or - s
      #s = float(sys.argv[1]) # Use the comman line arguement as specified in the R script
      #mut_map[sid] = freq_array[numpy.random.randint(0, len(freq_array))]*s
          
    md["selection_coeff"] = mut_map[sid]
  _ = tables.mutations.append(
          m.replace(metadata={"mutation_list": md_list})
  )

# check we didn't mess anything up
assert tables.mutations.num_rows == ots.num_mutations
print(f"The selection coefficients range from {min(mut_map.values()):0.2e}")
print(f"to {max(mut_map.values()):0.2e}.")



###### Correct the metadata and save ######

ts_metadata = tables.metadata
ts_metadata["SLiM"]["model_type"] = "WF"
tables.metadata = ts_metadata
ots = tables.tree_sequence()
ots.dump(output_path)
