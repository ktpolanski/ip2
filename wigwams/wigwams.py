import numpy as np
import multiprocessing as mp
import pandas as pd
import itertools as it
import csv
import scipy.stats as sps
import os
import time
import matplotlib.pyplot as plt
import seaborn as sns

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
	#computing the PCC is pretty simple and can be done in one shot
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
				#curr_pdf[j] = pvals_set[j] + np.log10(1 - np.power(10,(pvals_set[j+1]-pvals_set[j])))
				#numeric improvement. only matters for super ultra tiny cases, doesn't affect significance
				curr_pdf[j] = pvals_set[j] + np.log1p(-10**(pvals_set[j+1]-pvals_set[j])) / np.log(10)
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
	f.write(','.join(module[0]))
	
	#seed gene used, set size used and p-value obtained
	f.write('\t'+module[1])
	f.write('\t'+str(module[2]))
	f.write('\t'+str(module[3]))
	
	#gene membership, comma delimited
	f.write('\t'+','.join(module[4]))
	
	#new line
	f.write('\n')

def write_modules(fname, modules):
	'''
	A function that writes a whole np.array's worth of modules to a given file name.
	
	Input:
		* fname - the path where to store the modules
		* modules - the np.array of the modules to write
	'''
	
	writer = open(os.path.normcase(fname),'w')
	for i in range(modules.shape[0]):
		#sort the gene list, just in case
		modules[i,4].sort()
		write_module(writer,list(modules[i,:]))
	writer.close()

def read_modules(fname):
	'''
	A function that reads all of the modules from a given file and returns then as an np.array.
	
	Input:
		* fname - name of file to read the modules from
	'''
	
	#sponge up the file
	with open(fname, 'r') as f:
		reader = csv.reader(f,delimiter='\t')
		reader = list(reader)
	#now parse the bits and pieces. items 1 and 5 are actually lists
	for i in range(len(reader)):
		reader[i][0] = reader[i][0].split(',')
		reader[i][2] = int(reader[i][2])
		reader[i][3] = float(reader[i][3])
		reader[i][4] = reader[i][4].split(',')
	#that's about it. feed it out for other functions to use
	return np.asarray(reader)

def print_time(t1, t0, jobname):
	'''
	Helper function to print the time taken by procedures.
	'''
	m, s = divmod(round(t1-t0), 60)
	h, m = divmod(m, 60)
	holderstring = 'Took %d:%02d:%02d. '+jobname+' complete.'
	print(holderstring % (h, m, s))

def module_size_check(modules):
	'''
	Helper function that assesses the total module size and the number of unique genes.
	'''
	
	totsize = 0
	genes = []
	for i in range(modules.shape[0]):
		totsize += len(modules[i,4])
		genes = list(set(genes).union(set(modules[i,4])))
	#okay, we gots what we wants
	print('Number of modules: '+str(modules.shape[0]))
	print('Total module sizes: '+str(totsize))
	print('Number of unique genes: '+str(len(genes)))
	return(len(genes)/totsize)

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

def singlemining(seed, expr_df, deg_df, sets, alpha, corrnet, pvals, legacy):
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
		* legacy - toggle between most significant and biggest significant module output
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
		#note the [0][0] - np.where provides a bizarrely nested output and this digs out the actual value
		if legacy:
			ind = np.where(phold == np.min(phold))[0][0]
		else:
			#new idea to increase unique genes, pick biggest significant module
			#this will be the last p-value in the vector that clears the alpha
			ind = np.where(phold <= alpha)[0][-1]
		(overlap, pvalue) = singlemoduletest(corrgenes, sets[ind], comb, pvals[len(comb)-2][i], deg_df, degconds)
		#prepare module!
		condspan = degconds[np.asarray(comb)]
		out.append([list(condspan),seed,sets[ind],pvalue,overlap])
	#spit out the output
	return out

def _pool_init(expr_df_, deg_df_, sets_, alpha_, corrnet_, pvals_, legacy_):
	'''
	Initialising the "global" variables for the parallel module mining workers.
	'''
	global expr_df, deg_df, sets, alpha, corrnet, pvals, legacy
	expr_df = expr_df_
	deg_df  = deg_df_
	sets = sets_
	alpha = alpha_
	corrnet = corrnet_
	pvals = pvals_
	legacy = legacy_

def pool_singlemining(seed):
	'''
	The wrapper function for the parallelised mining, calling the actual mining function within.
	'''
	
	#take all the things that _pool_init made "global" for the pool running and call the actual function
	return singlemining(seed, expr_df, deg_df, sets, alpha, corrnet, pvals, legacy)

def mining(expr_df, deg_df, pool=1, sets=[50, 100, 150, 200, 250], alpha=0.05, corrnet=0.7, legacy=False, job='job'):
	'''
	The Wigwams analysis proper, mining for instances of dependent co-expression across multiple conditions.
	
	Input:
		* expr_df - column-multiindexed (treatment + time), standardised expression CSV file parsed into a PANDAS data frame.
		* deg_df - DEG status for each of the genes across each condition in expr_df, PANDAS data frame. If lacking DEG info, make a mock one with ones all across.
		* pool - how big to make the parallel pool for mining. Default: 1 (not parallel)
		* sets - a list of the top co-expressed gene set sizes for use in co-expression evaluation. Default: [50, 100, 150, 200, 250]
		* alpha - the significance threshold for use in statistical evaluation, Bonferroni corrected within the code. Default: 0.05
		* corrnet - if the top co-expressed genes stop being at least this correlated with the seed, sets stop being evaluated for that seed. Default: 0.7
		* legacy - True to run as in Polanski et al 2014, False to pick the longest statistically relevant module instead of the most statistically relevant. Default: False
		* job - job name for output naming purposes. Default: job
	'''
	
	#time the bugger
	t0=time.time()
	
	#set up output
	if not os.path.exists(job):
		os.makedirs(job)
	
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
	writer = open(os.path.normcase(job+'/raw_modules.tsv'),'w')
	if pool > 1:
		#this stuff is kind of confusing. explanation inbound!
		#_pool_init is the handle, pointing to the _pool_init function
		#this function will be called once for every pool creation
		#this will pass all the actual non-seed variables to the pools
		#the pools will subsequently be used to chew through all the seeds
		#(note the tuple after the pointer, those are the arguments that will be passed to the initialiser)
		p = mp.Pool(pool, _pool_init, (expr_df, deg_df, sets, alpha, corrnet, pvals, legacy))
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
			out = singlemining(seed, expr_df, deg_df, sets, alpha, corrnet, pvals, legacy)
			for module in out:
				write_module(writer, module)
				modcount +=1
				if np.mod(modcount,1000)==0:
					print('Found '+str(modcount)+' modules so far...')
	#okay, we're done, close the file
	print('Found '+str(modcount)+' modules in total')
	writer.close()
	
	#how long did we take?
	t1 = time.time()
	print_time(t1, t0, 'Module mining')

def find_modules_combination(modules, cat):
	'''
	Find indices of modules spanning a given condition combination.
	
	Input:
		* modules - module structure, np.array form
		* cat - list of conditions that the modules of interest are to span
	'''
	inds = []
	for i in range(modules.shape[0]):
		if modules[i,0]==cat:
			inds.append(i)
	inds = np.asarray(inds)
	return inds

def singlemerging(toggle, del_inds, raw_modules, mod_means, inds, overlap, meancorr, corrfilt, cat, i, expr_df):
	'''
	Helper function to reduce visual clutter of merging. Still kind of visually cluttered due to huge number of inputs though.
	'''
	
	#Flip the toggle to 0, if it hits 1 in the code then re-run the merging while.
	toggle = 0
	for j in range(i+1,len(inds)):
		#this is the module we'll try to merge into it
		#no point merging it if it's already merged into something
		if inds[j] in del_inds:
			continue
		#if the code is still here, we're eligible to check if module j can be merged into module i
		common_genes = list(set(raw_modules[inds[i],4]) & set(raw_modules[inds[j],4]))
		#will merging happen?
		mergetog = 0
		if float(len(common_genes))/float(min(len(raw_modules[inds[i],4]),len(raw_modules[inds[j],4]))) >= overlap:
			#sufficient membership overlap
			mergetog = 1
		else:
			#test the correlation of the means
			mergetog = 1
			for k in range(len(cat)):
				#this is the means structure, so just refer to modules as i and j
				if sps.pearsonr(mod_means[i][k], mod_means[j][k])[0] < meancorr:
					#correlation condition busted, merging not happening
					mergetog = 0
					break
		if mergetog == 0:
			#merging not happening
			continue
		#if we're here, merging is happening
		toggle = 1
		new_genes = list(set(raw_modules[inds[j],4])-set(common_genes))
		for gene in new_genes:
			#need to make sure that it's correlated enough
			add_gene = 1
			for k in range(len(cat)):
				if sps.pearsonr(mod_means[i][k], expr_df.loc[gene,cat[k]])[0] < corrfilt:
					#not correlated enough to the mean, abort the proceedings
					add_gene = 0
					break
			if add_gene == 1:
				#that means we're adding the gene. let's add the gene. gene. add gene.
				raw_modules[inds[i],4].append(gene)
		#post-merging stuff - note the module got merged, and recompute mean
		del_inds.append(inds[j])
		holder = []
		for cond in cat:
			holder.append(np.mean(expr_df.loc[raw_modules[inds[i],4],cond],axis=0))
		mod_means[i] = holder
	
	#all this stuff we need back as we may have changed it
	return (toggle, del_inds, raw_modules, mod_means)

def merging(expr_df, overlap=0.3, meancorr=0.9, corrfilt=0.8, condspan = None, which_file='raw', legacy=False, job='job'):
	'''
	The merging of redundant Wigwams modules spanning the same conditions.
	
	Input:
		* expr_df - PANDAS data frame of expression
		* overlap - what fraction of the smaller module size has to be made of the overlap to trigger merging. Default: 0.3
		* meancorr - how correlated the means of modules have to be to trigger merging independent of membership overlap. Default: 0.9
		* corrfilt - how correlated the genes have to be to the mean of the larger module to get moved over. Default: 0.8
		* condspan - if provided, only merge modules spanning exactly condspan conditions in total. Default: None (just merge everything)
		* which_file - which file to import and merge the modules from. Default: raw
		* legacy - True to run as in Polanski et al 2014, False to have an extra loop over the entire merging to assess the completeness of the procedure. Default: False
		* job - job name for output naming purposes. Default: job
	'''
	
	t0 = time.time()
	print('Commencing within-condition-span redundant module merging.')
	if condspan:
		print('Performing procedure for modules spanning '+str(condspan)+' conditions.')
	#let us read the modules first
	raw_modules = read_modules(os.path.normcase(job+'/'+which_file+'_modules.tsv'))
	#so, what categories will we need to look at?
	categories = np.unique(raw_modules[:,0])
	
	#super, super corner case - turns out that sometimes in merging, mergeable modules are constructed later on down
	#example - a large module spanning some conditions is made, it leaves some crumbs which don't quite fit
	#the crumbs are pieced together to form something that looks redundant to the original
	#as such, we need to loop merging until that sort of stuff is gone
	
	#dummy variable for deleted modules - full, to trigger loop
	del_inds = [1]
	while del_inds:
		#we're in the loop. del_inds means we have done some merging in the prior iteration
		#but now we blank it to actually evaluate stuff
		del_inds = []
		#alright. loop over the categories and get stuff done
		for cat in categories:
			#so, are we gonna even merge stuff in this category?
			if condspan:
				if len(cat) != condspan:
					continue
			#get the indices of the modules spanning that condition combination
			inds = find_modules_combination(raw_modules, cat)
			#print updates only if large lots to process
			if len(inds)>=100:
				print('Large module clump: '+str(len(inds))+' ('+', '.join(cat)+')')
			#sort the buggers on length
			lengths = []
			for i in range(len(inds)):
				lengths.append(len(raw_modules[inds[i],4]))
			indmask = np.asarray([i[0] for i in sorted(enumerate(lengths), reverse=True, key=lambda x:x[1])])
			inds = inds[indmask]
			#okay, now we have them largest to smallest
			#precompute module means
			mod_means = []
			for i in range(len(inds)):
				holder = []
				for cond in cat:
					holder.append(np.mean(expr_df.loc[raw_modules[inds[i],4],cond],axis=0))
				mod_means.append(holder)
			#commence merging!
			for i in range(len(inds)-1):
				#this is the module we're gonna try to merge the others with
				#no point merging if it's already merged into something
				if inds[i] in del_inds:
					continue
				#repeat merging until nothing merged
				toggle = 1
				while toggle==1:
					(toggle, del_inds, raw_modules, mod_means) = singlemerging(toggle, del_inds, raw_modules, mod_means, inds, overlap, meancorr, corrfilt, cat, i, expr_df)
		#wipe the modules to wipe if there are any
		if del_inds:
			raw_modules = np.delete(raw_modules,del_inds,axis=0)
			if legacy:
				#this just to emulate old times for old times sake
				del_inds=[]
	#Okay, so we're done with merging. Export!
	print('Merged down to '+str(raw_modules.shape[0])+' modules.')
	writer = write_modules(job+'/merged_modules.tsv',raw_modules)
	#how long did we take?
	t1 = time.time()
	print_time(t1, t0, 'Module merging')

def singlesweeping(raw_modules, del_inds, inds_sup, inds_sub, overlap):
	'''
	A helper function that executes individual instances of sweeping.
	'''
	
	for ind1 in inds_sup:
		#if it's swept already, we don't get to sweep with it
		if ind1 in del_inds:
			continue
		for ind2 in inds_sub:
			#same story. no need to re-sweep something that's already gone
			if ind2 in del_inds:
				continue
			#if we made it here, we can try to sweep
			#identify overlap
			genes = list(set(raw_modules[ind1,4]) & set(raw_modules[ind2,4]))
			#check if overlap big enough
			if float(len(genes)) / float(len(raw_modules[ind2,4])) > overlap:
				#module kill a smaller module. pow pow pow
				del_inds.append(ind2)
	#well, that's it
	return(del_inds)

def sweeping(overlap=0.5, condspan=None, which_file='merged', job='job'):
	'''
	The sweeping of redundant Wigwams modules across subsets of condition spans.
	
	Input:
		* overlap - what fraction of the gene membership of the module spanning fewer conditions that has to be made up of the overlap to trigger sweeping. Default: 0.5
		* condspan - if provided, only sweep from modules spanning exactly condspan conditions in total. Default: None (just sweep from everything)
		* which_file - which file to run sweeping on. Default: merged
		* job - job name for output naming purposes. Default: job
	'''
	
	t0 = time.time()
	print('Commencing inter-condition-span redundant module sweeping.')
	if condspan:
		print('Performing procedure for modules spanning '+str(condspan)+' conditions.')
	#let us read the modules first
	raw_modules = read_modules(os.path.normcase(job+'/'+which_file+'_modules.tsv'))
	#so, what categories will we need to look at?
	categories = np.unique(raw_modules[:,0])
	#dummy variable for deleted modules
	del_inds = []
	#this stuff above is all similar to what we had in merging.
	#now we need to sort the categories on length, to actually sweep with big ones
	lengths = []
	for i in range(len(categories)):
		lengths.append(len(categories[i]))
	indmask = np.asarray([i[0] for i in sorted(enumerate(lengths), reverse=True, key=lambda x:x[1])])
	categories = categories[indmask]
	#all right then. have the categories, can commence sweeping
	for i in range(len(categories)-1):
		#if we've hit length two, that's it
		if len(categories[i])==2:
			break
		#otherwise, if condspan, are we sweeping this combination?
		if condspan:
			if len(categories[i]) != condspan:
				continue
		#otherwise, time to get cooking
		inds_sup = find_modules_combination(raw_modules, categories[i])
		for j in range(i+1,len(categories)):
			#can we sweep this condition subset?
			if not set(categories[j]).issubset(set(categories[i])):
				#we can't
				continue
			#if we're here, we can
			inds_sub = find_modules_combination(raw_modules, categories[j])
			#let's run the sweeping in a separate function to be aesthetic about the indents
			del_inds = singlesweeping(raw_modules, del_inds, inds_sup, inds_sub, overlap)
	#sweeping complete. kick out swept modules
	if del_inds:
		swept_modules = np.delete(raw_modules,del_inds,axis=0)
	else:
		swept_modules = raw_modules
	print('Swept down to '+str(swept_modules.shape[0])+' modules.')
	writer = write_modules(job+'/swept_modules.tsv',swept_modules)
	#how long did we take?
	t1 = time.time()
	print_time(t1, t0, 'Module sweeping')

def thresholding(sizes, condspan=None, which_file='swept', job='job'):
	'''
	Filter the final Wigwams output to modules of desired sizes.
	
	Input:
		* sizes - an N-1-element list (where N is the number of conditions in your input), defining the desired sizes for modules from 2 to N in condition span.
		* condspan - if provided, only filter the modules spanning exactly condspan conditions. Default: None (just filter everything)
		* which_file - which file to import modules from and filter them on size. Default: swept
		* job - job name for output naming purposes. Default: job
	'''
	
	#not sure why I'm even timing this. completeness? but this is super fast.
	t0=time.time()
	print('Commencing filtering modules on size.')
	if condspan:
		print('Performing procedure for modules spanning '+str(condspan)+' conditions.')
	#let us read the modules first
	raw_modules = read_modules(os.path.normcase(job+'/'+which_file+'_modules.tsv'))
	#dummy variable for deleted modules
	del_inds = []
	#a-cracking we shall go
	for i in range(raw_modules.shape[0]):
		condlen = len(raw_modules[i,0])
		#condspan ahoy - can we filter this guy
		if condspan:
			if condlen != condspan:
				continue
		if len(raw_modules[i,4])<sizes[condlen-2]:
			del_inds.append(i)
	#well that was quick. kick out filtered modules
	if del_inds:
		filtered_modules = np.delete(raw_modules,del_inds,axis=0)
	else:
		filtered_modules = raw_modules
	print('Filtered down to '+str(filtered_modules.shape[0])+' modules.')
	writer = write_modules(job+'/filtered_modules.tsv',filtered_modules)
	#how long did we take?
	t1 = time.time()
	print_time(t1, t0, 'Module filtering')

def make_figure(modules, i, expr_df, cond_span, stand, job):
	'''
	A helper function that makes the figure(s) for a single module.
	'''
	
	for j in range(len(cond_span)):
		#so we're plotting. let's get plotting
		fig = plt.figure(figsize=(10,6))
		#how subplot'ty are we getting tonight
		#determine number of columns first
		if len(cond_span[j])>4:
			ncol = 3
		else:
			ncol = 2
		#number of rows follows naturally
		nrow = np.ceil(len(cond_span[j])/ncol)
		#prepare plotting area
		with sns.axes_style("ticks"):
			fig,axes = plt.subplots(int(nrow),int(ncol),sharey=True,figsize=(10,6))
		#go forth and plot
		for k,axs in enumerate(axes.flatten()):
			if not k < len(cond_span[j]):
				#we're not plotting. hide this
				axs.axis('off')
				continue
			#are we plotting XKCD pale red or XKCD midnight?
			coltog = '#03012d'
			if cond_span[j][k] in modules[i,0]:
				coltog = '#d9544d'
			#just plot things
			axs.set_title(cond_span[j][k])
			axs.margins(0.04)
			#extract time points
			tps = np.asarray(expr_df.loc[:,cond_span[j][k]].columns.get_level_values('Time').tolist())
			for ii in range(len(tps)):
				tps[ii] = float(tps[ii])
			for gene in modules[i,4]:
				axs.plot(tps, expr_df.loc[gene,cond_span[j][k]].values,color=coltog,linewidth=1)
			axs.set_xlim([min(tps),max(tps)])
		#label "outer" axes
		for axs in axes[-1,:]:
			axs.set_xlabel("Time")
		for axs in axes[:,0]:
			if stand:
				axs.set_ylabel("Standardised Expression")
			else:
				axs.set_ylabel("Expression")
		#whole plot complete
		fname = 'Module'+str(i+1).zfill(int(np.floor(np.log10(modules.shape[0]))+1))
		if len(cond_span)>1:
			fname = fname+'_plot'+str(j)
		plt.tight_layout()
		sns.despine()
		plt.savefig(os.path.normcase('plots-'+job+'/'+fname+'.svg'), bbox_inches='tight')
		#make new figure later
		plt.close("all")

def export_module(modules, i, annot, writer):
	'''
	A helper function that prints a single module to a file.
	'''
	
	for gene in modules[i,4]:
		#generic header stuff first
		writer.write('\n'+str(i+1)+'\t'+', '.join(modules[i,0])+'\t'+gene)
		#now then. if there's an annot, write more stuff
		if annot is not None:
			#make extra tab happen
			writer.write('\t')
			holder = annot[annot[:,0]==gene,1:]
			#we gotta check if we got a hit
			if holder.size > 0:
				#we got a hit b0ss. just in case make it be first hit
				writer.write('\t'.join(holder[0,:]))
			else:
				#no hit
				writer.write('UNMAPPED')

def export(expr_df, deg_df, which_file='filtered', annot_file=None, hyper=None, stand=True, job='job'):
	'''
	A function that takes the modules from the "module dump" format and makes a proper export.
	
	Input:
		* expr_df - PANDAS data frame of gene expression
		* deg_df - PANDAS data frame of DEG status, only used to extract proper condition order
		* which_file - which file to actually export. Default: filtered
		* annot - optional TSV with the expr_df identifiers in the first column, public identifiers in the second, and any additional information in third and onward. Default: None (off)
		* hyper - base URL to create Excel-friendly hyperlinks in the module export for each gene in a module. Provide public gene name placement with {gene} in the string. Default: None (off)
		* stand - Boolean informing whether the data was internally standardised or not, for axis labelling purposes. Default: True
		* job - job name for original output identification. Default: job
	'''
	
	t0=time.time()
	#start things off by sorting out the annot hyper situation
	if annot_file:
		#apparently the argparse thing opens the file by itself
		annot = csv.reader(annot_file,delimiter='\t')
		annot = np.asarray(list(annot))
		#capitalise the first column to match our names. leave second column be just in case
		for i in range(annot.shape[0]):
			annot[i,0] = annot[i,0].upper()
	else:
		annot = None
	if annot_file and hyper:
		annot2 = np.empty((annot.shape[0], annot.shape[1]+1), dtype=object)
		annot2[:,0] = annot[:,0]
		annot2[:,1] = annot[:,1]
		annot2[0,2] = 'Hyperlink'
		if annot.shape[1]>1:
			annot2[:,3:]=annot[:,2:]
		for i in range(annot2.shape[0]):
			annot2[i,2] = '=HYPERLINK("'+hyper.format(gene=annot2[i,1])+'", "LINK")'
		annot = annot2
	#okay, this should hopefully sort out the annotation situation. it's in and formatted
	
	#set up plot output folder
	if not os.path.exists('plots-'+job):
		os.makedirs('plots-'+job)
	#import appropriate modules
	modules = read_modules(os.path.normcase(job+'/'+which_file+'_modules.tsv'))
	#so, what are our conditions?
	conditions = np.asarray(deg_df.columns.tolist())
	#set up individual plot sizes. maximum of 9 plots per thing
	plot_counts = np.ceil(len(conditions)/9)
	plot_counts = np.zeros(int(plot_counts))+np.floor(len(conditions)/plot_counts)
	overhang = len(conditions)%len(plot_counts)
	for i in range(overhang):
		plot_counts[i] +=1
	cond_ranges = np.zeros(len(plot_counts)+1)
	cond_ranges[1:] = np.cumsum(plot_counts)
	cond_span=[]
	for i in range(len(plot_counts)):
		holder = conditions[int(cond_ranges[i]):int(cond_ranges[i+1])]
		cond_span.append(holder)
	#so, we know how the plots will be distributed if there's 10+ conditions.
	#this is fringe as hell, but you never know what people will feed into these things.
	
	#crack open a writer handle and prepare the header
	writer = open('exported_modules-'+job+'.tsv','w')
	writer.write('Module ID\tCondition Span\t')
	if annot_file:
		writer.write('\t'.join(annot[0,:]))
	else:
		writer.write('Gene Identifier')
	for i in range(modules.shape[0]):
		print('Exporting module '+str(i+1)+'...')
		#we're printing the i'th module
		make_figure(modules, i, expr_df, cond_span, stand, job)
		#figure status - done
		export_module(modules,i,annot,writer)
	#there we go. that's a wrap. 
	writer.close()
	#how long did we take?
	t1 = time.time()
	print_time(t1, t0, 'Module export')