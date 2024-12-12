## Modified from: https://tskit.dev/pyslim/docs/stable/vignette_coalescent_diversity.html


######## Python script that uses msprime and pyslim to generate a neutral burnin simulation ########

import pyslim
import msprime
import numpy
import sys

############## Command line arguments #################

Ne = float(sys.argv[1])
n_ind = float(sys.argv[2])
sequence_length = float(sys.argv[3])
r_msp = float(sys.argv[4])
mu_msp = float(sys.argv[5])
shape = float(sys.argv[6])
scale = float(sys.argv[7])
msprime_output_path = sys.argv[8]
mut_ratio = float(sys.argv[9])
DFE = sys.argv[10]                     # Can be "n" (normal) or "g" (gamma)
mean_alpha = float(sys.argv[11])
SD_alpha = float(sys.argv[12])
Set_ID = sys.argv[13]
sim = sys.argv[14]

###### Simulate a neutral ancestry using msprime ######

demog_model = msprime.Demography()
demog_model.add_population(initial_size=Ne) # Ne

ots = msprime.sim_ancestry(samples=n_ind, sequence_length=sequence_length, demography=demog_model, recombination_rate=r_msp)

# Add annotations for SLiM

ots = pyslim.annotate(ots, model_type="WF", tick=1, stage="late")

breakpoints = ots.breakpoints(as_array=True)
print("There are", ots.num_trees, "trees, associated with breakpoints", breakpoints)
print("There are", len(ots.samples()), "samples")

####### Add mutations #######

# Specifies that these mutations are going to be drawn using the infinite alleles model and are going to be assigned type "m2" in SLiM

mut_model = msprime.SLiMMutationModel(type=2)

ots = msprime.sim_mutations(ots, rate = mu_msp, model = mut_model, keep = True)

print(f"The tree sequence now has {ots.num_mutations} mutations, at "
      f"{ots.num_sites} distinct sites.")


####### Assign selection coefficients #########

tables = ots.tables
tables.mutations.clear()

# Decide the relative frequencies of neutral vs selected mutations


freq_deleterious = 10000
freq_beneficial = int(freq_deleterious*mut_ratio)
freq_neutral = 0

# Creae a list of 0s and 1s and -1s in frequencies of neutral, beneficial, and deleterious mutations
# This list will be used to generate a random 0 or 1 or -1 to decide whether the mutation becomes neutral or beneficial or deleterious

freq_array = [0]*freq_neutral + [1]*freq_beneficial + [-1]*freq_deleterious

mut_map = {}
for m in ots.mutations():
  md_list = m.metadata["mutation_list"]
  slim_ids = m.derived_state.split(",")
  assert len(slim_ids) == len(md_list)
  for sid, md in zip(slim_ids, md_list):
    
    if sid not in mut_map:
      
      if DFE == "g":
        
        # DFE is gamma
        mut_map[sid] = freq_array[numpy.random.randint(0, len(freq_array))]*(numpy.random.gamma(shape = shape, scale = scale)) 

      if DFE == "n":
        mut_map[sid] = numpy.random.normal(mean_alpha, SD_alpha)
        
          
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
ots.dump(f"{msprime_output_path}/{Set_ID}_sim{sim}_neutral_burnin.trees")
