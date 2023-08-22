## Modified from: https://tskit.dev/pyslim/docs/stable/vignette_coalescent_diversity.html


######## Python script that uses msprime and pyslim to generate a neutral burnin simulation ########

import pyslim
import msprime
import numpy

###### Simulate a neutral ancestry using msprime ######

demog_model = msprime.Demography()
demog_model.add_population(initial_size=10000)

ots = msprime.sim_ancestry(samples=1000, sequence_length = 10000000, demography=demog_model, recombination_rate=1e-8)

# Add annotations for SLiM

ots = pyslim.annotate(ots, model_type="WF", tick=1, stage="late")

####### Add mutations #######

# Specifies that these mutations are going to be drawn using the infinite alleles model and are going to be assigned type "m2" in SLiM

mut_model = msprime.SLiMMutationModel(type=2)

ots = msprime.sim_mutations(ots, rate = 1e-9, model = mut_model, keep = True)

print(f"The tree sequence now has {ots.num_mutations} mutations, at "
      f"{ots.num_sites} distinct sites.")


####### Assign selection coefficients #########

tables = ots.tables
tables.mutations.clear()

# Decide the relative frequencies of neutral vs selected mutations

freq_neutral = 20
freq_beneficial = 0
freq_deleterious = 0

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
        
          # DFE is gamma
          mut_map[sid] = freq_array[numpy.random.randint(0, len(freq_array))]*(min(numpy.random.gamma(shape = 0.2, scale = 0.25), 1.0)) # The max() function prevents negative fitnesses

          # DFE is just + or - s
          #s = 0.03
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
ots.dump("C:/Academics/Post-doc/Va_simulations/3_With_neutral_burnin_msprime/b_Interim_files/Msprime_outputs/neutral_burnin.trees")
