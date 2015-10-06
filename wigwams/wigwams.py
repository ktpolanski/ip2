import numpy as np
import multiprocessing as mp
import pandas as pd
import itertools as it

def row_pearson(mat, prof):
	'''
	A Pearson Correlation Coefficient function tuned for reasonable fastness!
	Also for what Wigwams wants, i.e. one profile vs. many profiles
	
	Input:
		* mat - NxM matrix of all gene expression profiles for a given condition
		* prof - the profile of the gene of interest (for that same condition obviously)
	'''
	
	#prepare data for doing things
	mat = mat-np.mean(mat,axis=1)[:,None]
	prof = prof-np.mean(prof)
	return np.sum(mat*prof,axis=1)/(np.sqrt(np.sum(mat*mat,axis=1))*np.sqrt(np.sum(prof*prof)))

def log_hypergeometric(x, G1, G2, U):
	"""
	The computation of the log-scale p-value of the pairwise hypergeometric test.
	Mirrored after Matlab function proposed by Meng et al. (2009).

	Input:
		* x - the size of the overlap of the two sets you're comparing
		* G1 - the size of the first set
		* G2 - the size of the second set
		* U - the size of the test universe (encompassing both sets and all other elements)
	"""
	
	#if any of the input arguments are zeros, then return log p value of 0
	if (x*G1*G2*U == 0):
		return 0
	
	#just gonna compile their logftable thing within the function and use that as needed
	#it would appear that the largest index of it they ever access is U+1
	logf = np.log(np.arange(1,U+2))
	logf = np.insert(logf,0,0)
	logf = np.cumsum(logf)
	
	#there we go. one logftable ready. can commence computation of actual teststuff
	minG = min(G1,G2)
	a = []
	for ind in range(1,(minG-x+2)):
		hold = x+ind-1
		a.append(-logf[hold]-logf[G1-hold]-logf[G2-hold]-logf[U+hold-G1-G2])
	a = np.asarray(a)
	amax = max(a)
	return logf[G1]+logf[U-G1]+logf[G2]+logf[U-G2]-logf[U]+amax+np.log(sum(np.exp(a-amax)))

def pval_sumtest(N, set, genes, pdf):
	"""
	A helper function for pvalue_preneration to avoid excessive loop nesting.
	"""
	
	#the components to sum up later
	comps = np.zeros(set-N+1)
	
	for i in range(N,set+1):
		comps[i-N] = pdf[i] + log_hypergeometric(N,i,set,genes)/np.log(10)
	
	#component status - acquired. time to add them up
	out = comps[0]
	for i in range(1,len(comps)):
		out = out + np.log10(1 + np.power(10,(comps[i]-out)))
	return out

def pvalue_pregeneration(set_sizes,genes,stress_count):
	"""
	The pre-generation of p-values for Wigwams to use when testing significance.
	
	Input:
		* set_sizes - the sizes of the sets that Wigwams will be evaluating for overlaps. Each gets its own p-values.
		* genes - how many total genes are present in the data (universe size)
		* stress_count - how many total stresses are present in the dataset (need to compute p-values for all subsets of 2:that-many)
	"""
	
	#dummy variable for output holding
	pvals_store = []
	
	if (stress_count == 2):
		#if we're in here, we're doing the original basic pairwise pre-generation
		#(this means there's no original pvals)
		pvals = []
		for i in set_sizes:
			curr_pvals = np.zeros(i+1)
			for j in range(0,i+1):
				curr_pvals[j] = log_hypergeometric(j,i,i,genes)/np.log(10)
			pvals_store.append(curr_pvals)
	else:
		#get the smaller overlap's p-values first. cheat a little, tell it there's less overlap
		pvals = pvalue_pregeneration(set_sizes,genes,stress_count-1)
		pvals_prev = pvals[-1]
		for i in range(0, len(set_sizes)):
			#decomposing into a pdf first
			set = set_sizes[i]
			pvals_set = pvals_prev[i]
			curr_pdf = np.zeros(set+1)
			for j in range(0,set):
				curr_pdf[j] = pvals_set[j] + np.log10(1 - np.power(10,(pvals_set[j+1]-pvals_set[j])))
			#the final element is the last test result
			curr_pdf[set] = pvals_set[set]
			#test proper time, everything is ready
			curr_pvals = np.zeros(set+1)
			for j in range(0,set+1):
				#this is the proper p-value, no adjusting needed
				curr_pvals[j] = pval_sumtest(j, set, genes, curr_pdf)
			#save the p-values
			pvals_store.append(curr_pvals)
	
	#save it to the output and pass it on
	pvals.append(pvals_store)
	return pvals

def write_module(f, module):
	'''
	A function that writes a single module to the file handle.
	
	Input:
		* f - file writer handle
		* module - the list of the single module to write to the file
	'''
	
	#condition span, comma delimited
	f.write(module[0][0])
	for i in range(1,len(module[0])):
		f.write(','+module[0][i])
	
	#seed gene used, set size used and p-value obtained
	#note the module[2][0] - no clue why, but the set is passed as a list
	f.write('\t'+module[1]+'\t'+str(module[2][0])+'\t'+str(module[3]))
	
	#gene membership
	f.write('\t'+str(module[4][0]))
	for i in range(1,len(module[4])):
		f.write(','+module[4][i])
	f.write('\n')

def singlemoduletest(corrgenes, singleset, comb, singlepvals, deg_df, degconds):
	'''
	A helper function to assess a single module
	'''
	#extract the overlap as the intersection of the top genes across the comb conditions
	overlap = set(corrgenes[range(singleset),comb[0]])
	for i in range(1,len(comb)):
		overlap = overlap & set(corrgenes[range(singleset),comb[i]])
	#clear out potential non-DE genes from overlap by using deg_df
	#has to be DE everywhere, that means the sum of the rows of the deg_df has to be len(comb)
	holder = deg_df.loc[overlap, degconds[np.asarray(comb)]]
	overlap = holder[np.sum(holder,axis=1)==len(comb)].index.tolist()
	#the seed gene invariably shows up in the overlap, account for that by shifting by 1
	#(p-values start at an overlap of 0, which is essentially what we're seeing in that case)
	return (overlap, singlepvals[len(overlap)-1])
	

def singlemining(seed, expr_df, deg_df, sets, alpha, corrnet, pvals):
	'''
	The function that actually does the mining for a single seed gene.
	
	Input:
		* seed - the name of the seed gene we're mining for
		* expr_df - parsed Pandas expression df
		* deg_df - parsed Pandas DEG df
		* sets - the sets to use in the testing
		* alpha - the significance to use in the testing
		* corrnet - the correlation threshold to cut off mining
		* pvals - the pre-generated p-values
	'''
	#So, where's our seed DEG in?
	conditions = deg_df.columns.tolist()
	mask = (deg_df.loc[seed,:]==1).tolist()
	#you need to turn the iterator into a list, and then that into an array
	degconds = np.asarray(list(it.compress(conditions, mask)))
	#holder of gene names
	genes = expr_df.index.tolist()
	#pre-compute correlation matrix for relevant stresses, and strip it down
	corrmat = np.zeros((len(sets),len(degconds)))
	corrgenes = np.ndarray((np.max(sets),len(degconds)),dtype=object)
	for cond in range(len(degconds)):
		#compute the proper correlation matrix
		corrhold = row_pearson(expr_df.loc[:,degconds[cond]].values,expr_df.loc[seed,degconds[cond]].values)
		#now sort the bugger, also getting out indices of things
		corrinds = np.asarray([i[0] for i in sorted(enumerate(corrhold), reverse=True, key=lambda x:x[1])])
		corrhold = np.asarray(sorted(corrhold, reverse=True))
		#we only need the top max(sets) number of genes for later
		for i in range(0,np.max(sets)):
			corrgenes[i,cond]=genes[corrinds[i]]
		#in fact, we only need values of the correlation at the set count points
		#storage optimisation!
		for i in range(len(sets)):
			#note the -1, python indexing shift
			corrmat[i,cond]=corrhold[sets[i]-1]
	
	#prepare combinations - these come as indices for ease at corrmat pointing.
	#for module export, don't forget to call degconds of the appropriate indices.
	comblist = []
	for combcount in range(2,len(degconds)+1):
		comblist = comblist + list(it.combinations(range(len(degconds)),combcount))
	
	#okay, all the things are prepared now. time for testing proper
	out=[]
	for comb in comblist:
		phold = np.zeros(len(sets))
		for i in range(len(sets)):
			#check the net - if any of the top (set) genes across the conditions are below it, abort
			if any(corrmat[i,comb]<corrnet):
				break
			#otherwise, propose a module - notice the pvals indexing.
			#The first one shifted by 2 - as 2 sets is pvals[0]
			#The second one is just which set we're running this with
			(overlap, phold[i]) = singlemoduletest(corrgenes, sets[i], comb, pvals[len(comb)-2][i], deg_df, degconds)
		#okay, so we ran all the sets. check if there's anything relevant
		if all(phold > alpha):
			#sometimes there isn't
			continue
		#if we're here, then there is something worth it. let's dig it out!
		ind = np.where(phold == np.min(phold))
		(overlap, pvalue) = singlemoduletest(corrgenes, sets[ind], comb, pvals[len(comb)-2][i], deg_df, degconds)
		#prepare module!
		condspan = degconds[np.asarray(comb)]
		out.append([list(condspan),seed,sets[ind],pvalue,overlap])
	#spit out the output
	return out

def _pool_init(expr_df_, deg_df_, sets_, alpha_, corrnet_, pvals_):
	'''
	Initialising the "global" variables for the parallel module mining workers.
	'''
	global expr_df, deg_df, sets, alpha, corrnet, pvals
	expr_df = expr_df_
	deg_df  = deg_df_
	sets = sets_
	alpha = alpha_
	corrnet = corrnet_
	pvals = pvals_

def pool_singlemining(seed):
	'''
	The wrapper function for the parallelised mining, calling the actual mining function within.
	'''
	#take all the things that _pool_init made "global" for the pool running and call the actual function
	return singlemining(seed, expr_df, deg_df, sets, alpha, corrnet, pvals)

def mining(expr_df, deg_df, pool=1, sets=[50, 100, 150, 200, 250], alpha=0.05, corrnet=0.7):
	'''
	The Wigwams analysis proper, mining for instances of dependent co-expression across multiple conditions.
	
	Input:
		* expr_df - column-multiindexed (treatment + time), standardised expression CSV file parsed into a PANDAS data frame.
		* deg_df - DEG status for each of the genes across each condition in expr_df, PANDAS data frame. If lacking DEG info, make a mock one with ones all across.
		* pool - how big to make the parallel pool for mining. Default: 1 (not parallel)
		* sets - a list of the top co-expressed gene set sizes for use in co-expression evaluation. Default: [50, 100, 150, 200, 250]
		* alpha - the significance threshold for use in statistical evaluation, Bonferroni corrected within the code. Default: 0.05
		* corrnet - if the top co-expressed genes stop being at least this correlated with the seed, sets stop being evaluated for that seed. Default: 0.7
	'''
	
	#sort sets just in case, and array them
	sets = np.asarray(sorted(sets))
	
	#pre-generation of seeds
	seeds = deg_df[deg_df.sum(axis=1)>1].index.tolist()
	
	#pre-generate p-values
	print('Pre-generating p-values...')
	pvals = pvalue_pregeneration(sets, len(deg_df.index), len(deg_df.columns))
	
	#fix up the alpha - make it be log10 scale and Bonferroni'd
	tcount = np.sum(deg_df.values, axis=1)
	tcount = np.power(2,tcount) - tcount - 1
	alpha = np.log10(alpha / np.sum(tcount))
	
	#mining proper
	print('Commencing module mining...')
	modcount = 0
	writer = open('raw_modules.txt','w')
	if pool > 1:
		#this stuff is kind of confusing. explanation inbound!
		#_pool_init is the handle, pointing to the _pool_init function
		#this function will be called once for every pool creation
		#the pools will subsequently be used to chew through all the seeds
		#(note the tuple after the pointer, those are the arguments that will be passed to the initialiser)
		p = mp.Pool(pool, _pool_init, (expr_df, deg_df, sets, alpha, corrnet, pvals))
		#imap_unordered makes the output unordered, but provides it fast
		#this is an iterator, so you for over it to process the things as they come out
		for out in p.imap_unordered(pool_singlemining, seeds):
			for module in out:
				write_module(writer, module)
				modcount += 1
				if np.mod(modcount,1000)==0:
					print('Found '+str(modcount)+' modules so far...')
	else:
		for seed in seeds:
			out = singlemining(seed, expr_df, deg_df, sets, alpha, corrnet, pvals)
			for module in out:
				write_module(writer, module)
				modcount +=1
				if np.mod(modcount,1000)==0:
					print('Found '+str(modcount)+' modules so far...')
	#okay, we're done, close the file
	print('Found '+str(modcount)+' modules in total')
	writer.close()