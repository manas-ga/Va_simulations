# Read the output file created by SliM
f_simout = open("SLiMeuro_out.txt").read()

# Create an empty file to write the mutations in
f_mutation = open("mutations.txt", "w")

dict_1 = {}
dict_2 = {}

for i, line in enumerate(open("SLiMeuro_out.txt")):
    dict_1[line.rstrip()] = i
    dict_2[i] = line

# Identify start and end points of the section showinng mutations
start = dict_1["Mutations:"]
end = dict_1["Individuals:"]

#print(start)
#print(end)

# Write headers
f_mutation.write("Temp_ID Permanent_ID Position s h Subpop Tick Number\n")
for i in range(start + 1,end):
   f_mutation.write(dict_2[i])

f_mutation.close()



        


