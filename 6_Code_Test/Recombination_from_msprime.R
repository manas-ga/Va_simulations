

######################################################################################################
######################################################################################################
# Function to simulate the population
######################################################################################################
######################################################################################################

	# Nflies=100						#number of diploid flies. Assumed to be half males and half females
	# NGen=1							#number of generations
	# HapLength=10000					#number of BASES in the region
	# RecombRate=3e-8					#recomb rate per base per generation for base populations
	# ExpectedRecombs=HapLength*RecombRate	#expected recombinations for the forward simulation, or enter your own expected value if you want something higher
	# MinRecombsPerGamete=1				# or force a minimum number per haplotype
	# MutationRate=1.5e-9				#mutation rate per base per gen in MSprime
	# Ne=1e6							#Ne for MSprime
	# Fstem<-NA							#stem for the output file names. Can be left as NA
	# randomise_state_labels=TRUE		#Set False if derived mutations are coded 1
	# sda=0.1							#standard deviation of additive effect
	# sdwe=0.1							#standard deviation of environmental effect on fitness



simulate_base_population<-function(Nflies=10,HapLength=10000,RecombRate=3e-8, MutationRate=1.5e-9,Ne=1e6,Fstem=NA, randomise_state_labels=TRUE, msp_path=""){

	#make a file name if required
	if(is.na(Fstem)){Fstem<-gsub(" ","_",paste(date(),paste(sample(0:9,10,replace=TRUE),collapse=""),collapse=""))}
	
	print("Simulating base Population")
	
	#simulate the starting population using msprime
	sim_filename<-paste(Fstem,".msp",sep="")
	paste(paste0(msp_path, "msp "), "simulate --length ",HapLength," --recombination-rate ",RecombRate," --mutation-rate ",MutationRate," --effective-population-size ",Ne," ",2*Nflies," ",sim_filename,sep="")->mspString
	print(mspString)
	system(mspString)
	vars_filename<-paste(Fstem,".msp.vars",sep="")
	system(paste(paste0(msp_path, "msp "), "variants ",sim_filename," > ",vars_filename,sep=""))
	haps_filename<-paste(Fstem,".msp.haps",sep="")
	system(paste(paste0(msp_path, "msp "), "haplotypes ", sim_filename," > ",haps_filename,sep=""))
	#read in our starting haplotypes
	round(read.table(vars_filename,sep="\t")[,1],0)->vars
	as.matrix(read.fwf(haps_filename,widths=rep(1,length(vars))))->starting_haps
	#get inter-varient relative genetic distances # We could multiply this by the rate, if we had a map of rate variation that was used to simulate the base population
	
	#randomise the 1/0 labelling here
	if(randomise_state_labels){
	  sample(0:1,ncol(starting_haps),replace=TRUE)->state_modifier
	  t(apply(starting_haps,1,function(x){(x+state_modifier)%%2}))->starting_haps
	}

	
	return(list(Parental=starting_haps,SNPs=vars,mspFile=Fstem))
	
	}
	
simulate_forward_population<-function(NGen=2,fitness_effects,starting_haps,vars,HapLength=10000,RecombRate=3e-8,ExpectedRecombs=NA, AtleastOneRecomb=FALSE, Fstem=NA, sdwe=0.1){


	##################################################################	
	#parental haplotypes in starting_haps matrix
	#SNP locations in vars vector
	##################################################################
	
	Nflies<-nrow(starting_haps)/2
	if(is.na(ExpectedRecombs)){ExpectedRecombs<-HapLength*RecombRate}
	
	#find the location of our varients, and shift the overlapping ones to become neighbours (**BAD HACK**)
	while(max(table(vars))>1){vars[duplicated(vars)]<-vars[duplicated(vars)]+1}
	length(vars)->nVars	

	vars[2:nVars]-vars[1:(nVars-1)]->var_distances

	#set up the matrices for haplotypes. I hope that setting them up here leads to less memory allocation in the loop
	new_generation<-old_generation<-starting_haps
	parental_genotypes<-matrix(0,ncol=nVars,nrow=Nflies)
	
	#set up to record where the recombinations events happen
	recomb_events<-0
	recomb_locations<-vector()

	print("Simulating Forwards")	


	##################################################################	
	# Preliminary generation in which each wild mother and father contribute 1 male and 1 female to the starting population, without selection
	##################################################################	

	for (o in 1:Nflies){

		o_indices<-c(2*o-1,2*o)  #where the offspring haplotypes will be located
		
		#Maternal contribution. I assume the first half of the flies are female, and the second half of the flies are male 
		#Female offspring N will have motehr 1 and father 1
		mother<-((o-1)%%(Nflies/2))+1
		#locate her haplotypes
		mat_indices<-sample(c(2*mother-1,2*mother) )#where the mothers haplotypes are located, but in a random order
		#select the number of recombination events

		if(AtleastOneRecomb){
			mat_recomb_events<-actuar::rztpois(1,ExpectedRecombs)
		}else{
			mat_recomb_events<-rpois(1,ExpectedRecombs)
		}	
		#if there is recombination

		if(mat_recomb_events>0){
			sort(sample(1:(nVars-1),mat_recomb_events,prob=var_distances))->where #which SNPs do recombination occur to the immediate right of - allowing for distances between them
			recomb_events<-recomb_events+mat_recomb_events	#running total of the number
			recomb_locations<-c(recomb_locations,where) #running list of the locations. Bad memory usage, but the numebrs are small!
			unlist(mapply(rep,rep(c(TRUE,FALSE), length=mat_recomb_events+1),c(where[1],where[-1]-where[-mat_recomb_events],nVars-where[mat_recomb_events])))->haplo_pattern # excersize for the reader! repeat TRUE for all vars up to first recomb, then false to the next, and so on....
			#remembering that the order of the maternal haplotype indices has been randomised, so which haplotype starts is random
			new_generation[o_indices[1],haplo_pattern]<-old_generation[mat_indices[1],haplo_pattern]
			new_generation[o_indices[1],!haplo_pattern]<-old_generation[mat_indices[2],!haplo_pattern]
			
		}else{
		#else if no recombination, pick one of her haplotypes at random
			new_generation[o_indices[1],]<-old_generation[mat_indices[1],] #note that the order of the mat indices is already randomised
		}
		
		#Paternal contribution. I assume the first half of the flies are female, and the second half of the flies are male 
		#Male offspring N will have mother 1 and father 1
		father<-(Nflies/2)+(((o-1)%%(Nflies/2))+1)
		#locate his haplotypes
		pat_indices<-sample(c(2*father-1,2*father) )#where the fathers haplotypes are located,but in a random order
		#select the number of recombination events
		if(AtleastOneRecomb){
			pat_recomb_events<-actuar::rztpois(1,ExpectedRecombs)
		}else{
			pat_recomb_events<-rpois(1,ExpectedRecombs)
		}
		#if there is recombination
		if(pat_recomb_events>0){
			sort(sample(1:(nVars-1),pat_recomb_events,prob=var_distances))->where #which SNPs do recombination occur to the immediate right of - allowing for distances between them
			recomb_events<-recomb_events+pat_recomb_events	#running total of the number
			recomb_locations<-c(recomb_locations,where) #running list of the locations. Bad memory usage, but the numebrs are small!
			unlist(mapply(rep,rep(c(TRUE,FALSE), length=pat_recomb_events+1),c(where[1],where[-1]-where[-pat_recomb_events],nVars-where[pat_recomb_events])))->haplo_pattern # excersze for the reader! repeat TRUE for all vars up to first recomb, then false to the next, and so on....
			#remembering that the order of the paternal haplotype indices has been randomised, so which haplotype starts is random
			new_generation[o_indices[2],haplo_pattern]<-old_generation[pat_indices[1],haplo_pattern]
			new_generation[o_indices[2],!haplo_pattern]<-old_generation[pat_indices[2],!haplo_pattern]
			
		}else{
		#else if no recombination, pick one of his haplotypes at random
			new_generation[o_indices[2],]<-old_generation[pat_indices[1],] #note that the order of the pat indices is already randomised
		}
		
	} #end of (For each fly)
	
	FirstGeneration<-old_generation<-new_generation

	##################################################################	
	# Then simulate forwards for N generations, with selection
	##################################################################	
		
	#Additive genetic variation in fitness
	 # effect of each locus on fitness  # could associate them with frequencies here if we want
	
	#for each generation
	for(g in 1:NGen){
		
		if((g%%10)==0){print(paste("Generation",g))}
		
		#model the environmental fitnesses
		we<-rnorm(Nflies,0,sdwe) # fitness due to environment# 
		
		#parental_genotypes is a matrix with loci in columns and individuals in rows. 
		#created here #### It takes the number of reference alleles (which are arbitrarily defined)
		parental_genotypes<-(old_generation[seq(1,2*Nflies,by=2),]+old_generation[seq(2,2*Nflies,by=2),])/2
		parental_fitnesses<-exp(parental_genotypes%*%fitness_effects+we) 
		
		#for each new offspring fly
		for (o in 1:Nflies){
			o_indices<-c(2*o-1,2*o)  #where the offspring haplotypes will be located
			
		#Maternal contribution. I assume the first half of the flies are female, and the second half of the flies are male 
			#select a random mother in proportion to her relative fitness
			mother<-sample(1:(Nflies/2),size=1,prob=parental_fitnesses[1:(Nflies/2)])
			#locate her haplotypes
			mat_indices<-sample(c(2*mother-1,2*mother) )#where the mothers haplotypes are located,but in a random order
			#select the number of recombination events

			if(AtleastOneRecomb){
				mat_recomb_events<-actuar::rztpois(1,ExpectedRecombs)
			}else{
				mat_recomb_events<-rpois(1,ExpectedRecombs)
			}

			#if there is recombination
			if(mat_recomb_events>0){
				sort(sample(1:(nVars-1),mat_recomb_events,prob=var_distances))->where #which SNPs do recombination occur to the immediate right of - allowing for distances between them
				recomb_events<-recomb_events+mat_recomb_events	#running total of the number
				recomb_locations<-c(recomb_locations,where) #running list of the locations. Bad memory usage, but the numebrs are small!
				unlist(mapply(rep,rep(c(TRUE,FALSE), length=mat_recomb_events+1),c(where[1],where[-1]-where[-mat_recomb_events],nVars-where[mat_recomb_events])))->haplo_pattern # excersize for the reader! repeat TRUE for all vars up to first recomb, then false to the next, and so on....
				#remembering that the order of the maternal haplotype indices has been randomised, so which haplotype starts is random
				new_generation[o_indices[1],haplo_pattern]<-old_generation[mat_indices[1],haplo_pattern]
				new_generation[o_indices[1],!haplo_pattern]<-old_generation[mat_indices[2],!haplo_pattern]
				
			}else{
			#else if no recombination, pick one of her haplotypes at random
				new_generation[o_indices[1],]<-old_generation[mat_indices[1],] #note that the order of the mat indices is already randomised
			}
			
		#Paternal contribution. I assume the first half of the flies are female, and the second half of the flies are male 
			#select a random father in proportion to his relative fitness
			father<-sample(((Nflies/2)+1):Nflies,size=1,prob=parental_fitnesses[((Nflies/2)+1):Nflies])
			#locate his haplotypes
			pat_indices<-sample(c(2*father-1,2*father) )#where the fathers haplotypes are located,but in a random order
			#select the number of recombination events
			if(AtleastOneRecomb){
				pat_recomb_events<-actuar::rztpois(1,ExpectedRecombs)
			}else{
				pat_recomb_events<-rpois(1,ExpectedRecombs)
			}
			#if there is recombination
			if(pat_recomb_events>0){
				sort(sample(1:(nVars-1),pat_recomb_events,prob=var_distances))->where #which SNPs do recombination occur to the immediate right of - allowing for distances between them
				recomb_events<-recomb_events+pat_recomb_events	#running total of the number
				recomb_locations<-c(recomb_locations,where) #running list of the locations. Bad memory usage, but the numebrs are small!
				unlist(mapply(rep,rep(c(TRUE,FALSE), length=pat_recomb_events+1),c(where[1],where[-1]-where[-pat_recomb_events],nVars-where[pat_recomb_events])))->haplo_pattern # excersze for the reader! repeat TRUE for all vars up to first recomb, then false to the next, and so on....
				#remembering that the order of the paternal haplotype indices has been randomised, so which haplotype starts is random
				new_generation[o_indices[2],haplo_pattern]<-old_generation[pat_indices[1],haplo_pattern]
				new_generation[o_indices[2],!haplo_pattern]<-old_generation[pat_indices[2],!haplo_pattern]
				
			}else{
			#else if no recombination, pick one of his haplotypes at random
				new_generation[o_indices[2],]<-old_generation[pat_indices[1],] #note that the order of the pat indices is already randomised
			}
			
		} #end of (For each fly)
		
		old_generation<-new_generation
		
	}	#end of (For each generation)
	
	return(list(Generations=g,SNPs=nVars,Parental=starting_haps,FirstGeneration=FirstGeneration,LastGeneration=new_generation,Fstem=Fstem,TotalRecombinationEvents=recomb_events,AllRecombinationLocations=recomb_locations))
	
}



