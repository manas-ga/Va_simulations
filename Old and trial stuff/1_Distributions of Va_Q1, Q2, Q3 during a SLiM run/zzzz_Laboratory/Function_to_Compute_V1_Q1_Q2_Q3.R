#############################################################################
############ Details of the function ########################################
#############################################################################

## Input variables: 
# 1. Path of the file storing the csv file containing the c matrix for genomes (path_genomes)
# 2. Path of the text file storing the information on mutations (path_mutations)

## What does the function do?
# Reads the c matrix of genomes generated by the python code  (i.e. sampled individuals)
# Calculates the c matrix for individuals of the sample 
# Calculates the L matrix 
# Deletes loci that are not segregating in the sample and calculates the R matrix for the sample
# Reads the list of mutations generated by the python code and calculates fitness for every individual
# Fits a linear regression for Fitness~C
# Calculates Va, Q1, Q2, and Q3

## Output
# The function returns va, q1, q2, and q3


#specify the path of the file containing the genomes and mutations (outputs from SLim and Python)
#setwd("C:/Academics/Post-doc/SLiM/Trial stuff/Trial3 (alpha for non-additivity)")

########################################
######## Define the function ###########
########################################

#compute_va_q1_q2_q3 = function(path_genomes, path_mutations){

########################################
### Compute C, L and R matrices ########
########################################

path_mutations = "b_Interim_files/Mutations/mutations_176000.txt"
path_genomes = "b_Interim_files/C_Matrices/c_matrix_176000.csv"

# Read genomes
c_genome = read.csv(path_genomes, header=F)

#Convert genome data into individual data (rows are individuals and columns are allele counts {0,1 or 2} at various sites). Note that genome 1 and genome 2 are from individual 1, and so on

n_individuals = nrow(c_genome)/2
n_sites = ncol(c_genome)

c_ind = c_genome[seq(1, 2*n_individuals, 2),] + c_genome[seq(2, 2*n_individuals, 2),]

############## !!!!!!!!!!Test!!!!!!!!!!!!
## Removig two individuals with REALLY low fitnesses
#c_ind = c_ind[-c(10, 52),]



#### Calculate the matrix of second mixed moments of c, ie. the L matrix

L = cov(c_ind)

# Check for 0s on the diagonal of L
#which(diag(L)==0)

####### Notice that L has some 0s on the diagonal
# This is happening because when I sample individuals in SLiM from a population of 10000, the sample does not contain some low frequency mutations, ie. some loci are not segregating in the sample, but are in the original population

# identify the loci that are missing in the sample

missing_loci = which(diag(L)==0)
retained_loci = which(diag(L)!=0)

# Trim L to contain only the retained loci

L_ret = L[retained_loci,retained_loci]

# Trim the c matrix to contain only the retained loci, i.e. the loci that are segregating in the sample

c_ind_ret = c_ind[,retained_loci]

# Correlation matrix

R = cov2cor(L_ret)



##########################################
### Calculate fitnesses of individuals ###
##########################################

## Read the list of mutations output generated by SLiM and subsequently cropped out as a separate file by Python

mutations = read.table(path_mutations, header=T, sep = " ")


#mutations$Frequency = mutations$Number/20000
#mutations$Diversity = 2*(mutations$Frequency)*(1-mutations$Frequency)


# Sort mutations based on the Temp_ID, so that the order of loci matches the order in the C and L matrices

mutations = mutations[order(mutations$Temp_ID),]

# Trim mutations to contain only the retained loci, i.e. the loci that are segregating in the sample

mutations_ret = mutations[retained_loci,]

# Calculate the number of loci that are retained

n_sites_ret = ncol(c_ind_ret)

# Create an empty vector of individual fitnesses with each element initialized to 1 
#This will have to be 0 when using an additive model for fitness

w = rep(1, nrow(c_ind_ret))

# Calculate the fitness of each individual (i.e. populate w) by multiplying the contributions of each locus
# The product will have to be replaced by a sum when things are additive
# k loops over individuals, and m loops over sites

for (k in 1:n_individuals){
  # loop through loci, summing fitness contributions
  for (m in 1:n_sites_ret){
    # Calculate the fitness of each individual using a multiplicative model for fitness
    w[k] = w[k] * (1 + ((2-c_ind_ret[k,m])*mutations_ret$h[m] + (c_ind_ret[k,m]-1)/2)*c_ind_ret[k,m]*mutations_ret$s[m]) 
  
    }
}


###################################
### Perform linear regression #####
###################################

#Calculate relative
Fitness = w/mean(w)

#TEST BY REMOVING INDIVIDUALS WITH EXCEPTIONALLY LOW FITNESSES
#fit_1 <- lm(Fitness[-c(10, 52)] ~ ., data = c_ind_ret[-c(10, 52),])

fit_1 <- lm(Fitness ~ ., data = c_ind_ret)


######### The section on reordering is not required for this function #########################################

# Since n_sites > n_individuals, all coefficients (ie, all alphas) cannot be estimated
# There is aliasing
# To investigate if the aliasing affects calculations of V_a, reorder the data frame used for the linear model

#vect_reordered = sample(1:n_sites_ret, n_sites_ret, replace = F)
#vect_reordered1 = sample(1:n_sites_ret, n_sites_ret, replace = F)


#fit_reordered = lm(Fitness ~ ., data = c_ind_ret[,vect_reordered])
#fit_reordered1 = lm(Fitness ~ ., data = c_ind_ret[,vect_reordered1])

#fit_2 <- lm(w ~ ., data = c_ind)
# Notice that fit_1 and fit_2 produce subtly different outputs


#summary(fit_1)
#summary(fit_2)
#summary(fit_reordered)

#################################################################################

# Remove the first coefficient (i.e. the intercept), and store the alphas

list_alpha = coef(fit_1)[-1] 
#list_alpha_reordered = coef(fit_reordered)[-1] 
#list_alpha_reordered1 = coef(fit_reordered1)[-1] 


# Since n_individuals < n_loci, many coefficients cannot be estimated, and list_alpha, etc have NAs

#Remove NAa

list_alpha[is.na(list_alpha)]<-0
#list_alpha_reordered[is.na(list_alpha_reordered)]<-0
#list_alpha_reordered1[is.na(list_alpha_reordered1)]<-0


##### Some sanity checks #####

#list_freq = mutations$Frequency
#list_freq1 = 0.5*colMeans(c_ind)
#plot(list_freq, list_freq1)
#plot(list_s, list_alpha)
#plot(log(mutations$Diversity), log(mutations$s*mutations$s))
#lm(log(list_alpha*list_alpha)~log(2*list_freq*(1-list_freq)))
#lm(log(mutations$s*mutations$s)~log(mutations$Diversity))


################################################################
######### Quantities to be estimated ###########################
################################################################

##### Buffalo and Coop (2019) treat selected and neutral loci separately
## In order to test some of the assumptions made by B&C, we need to compute quantities Q1, Q2, and Q3
# This requires summing over selected and neutral loci independently

## IDs of neutral and selected mutations

list_neutral = which(mutations_ret$s==0)
list_selected = which(mutations_ret$s!=0)

#list_neutral_reordered = which(mutations_ret[vect_reordered,]$s==0)
#list_selected_reordered = which(mutations_ret[vect_reordered,]$s==0)


###########
### va ####
###########

# Calculate Va

va = t(list_alpha)%*%L_ret%*%list_alpha
#va_reordered = t(list_alpha_reordered)%*%cov(c_ind_ret[,vect_reordered])%*%list_alpha_reordered
#va_reordered1 = t(list_alpha_reordered1)%*%cov(c_ind_ret[,vect_reordered1])%*%list_alpha_reordered1

# Does Va remain invariant to aliasing?
#va
#va_reordered
#va_reordered1

###########
### Q1 ####
###########

# Calculate frequencies at the retained loci

Frequency = colMeans(c_ind_ret)/2

q1 = cov(2*0.5*Frequency*(1-Frequency), list_alpha*list_alpha)

#q1_reordered = cov(2*0.5*Frequency[vect_reordered]*(1-Frequency[vect_reordered]), list_alpha_reordered*list_alpha_reordered)
#q1_reordered1 = cov(2*0.5*Frequency[vect_reordered1]*(1-Frequency[vect_reordered1]), list_alpha_reordered1*list_alpha_reordered1)

#General expression (to check if reordering loci affects Q1). Running the following command gives a different output every time
#q1
#q1_reordered
#q1_reordered1

# Q1 is not invariant to aliasing!

###########
### Q2 ####
###########

# Loop over pairs of selected loci using j and k. The sum function loops over neutral loci

# Initialise q2 to 0, and then recursively sum over pairs of values of j and k

q2 = 0

for (j in 1:(length(list_selected) - 1)){
  for (k in (j+1):length(list_selected)){
    q2 = q2 + list_alpha[list_selected[j]]*list_alpha[list_selected[k]]*((2*Frequency[list_selected[j]]*(1 - Frequency[list_selected[j]])*2*Frequency[list_selected[k]]*(1 - Frequency[list_selected[k]]))^0.5)*sum(R[j, list_neutral]*R[k, list_neutral])
    
  }
  #print(q2)
}

#Q2 


###########
### Q3 ####
###########

# Q3 is a covariance. Calculating the two terms (q3_a and q3_b) separately

q3_a = rowSums(R[list_selected, list_neutral]*R[list_selected, list_neutral])
q3_b = list_alpha[list_selected]*((2*Frequency[list_selected]*(1 - Frequency[list_selected]))^0.5)

q3 = cov(q3_a, q3_b)

# Using reordered loci to check if q3 is invariant to aliasing

#q3_a_reordered = rowSums(R[list_selected_reordered, list_neutral_reordered]*R[list_selected_reordered, list_neutral_reordered])
#q3_b_reordered = list_alpha_reordered[list_selected_reordered]*((2*Frequency[vect_reordered][list_selected_reordered]*(1 - Frequency[vect_reordered][list_selected_reordered]))^0.5)

#q3_reordered = cov(q3_a_reordered, q3_b_reordered)

#q3
#q3_reordered

##############################
######### Output #############
##############################

#return(c(va, q1, q2, q3))

#}

