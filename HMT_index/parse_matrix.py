import argparse
import numpy as np
import pandas as pd
import scipy.stats as sps
import sys

def overlapcheck(span1, span2):
	#helper function to assess whether span1 and span2 overlap
	if span1[0] >= span2[0] and span1[0] <= span2[1]:
		return True
	if span1[1] >= span2[0] and span1[1] <= span2[1]:
		return True
	return False

def hommel(pvals):
	#helper function to combine multiple p-values into a single p-value
	#not using Fischer method because that one just makes the p-values tiny
	#when I want a sort of "mean/median" thing
	#taken from Hommel (1983) formula reproduced in Vovk (2012)
	
	#pre-processing - is our thing a single value? if so, just return it
	if isinstance(pvals, float):
		return pvals
	#well, if it's not a single value, we can get going. turn to np.array and sort
	pvals = np.sort(np.asarray(pvals))
	#actually computing the thing
	holder = np.arange(len(pvals))+1
	pvals = (pvals * len(pvals)) / holder
	return (np.sum(1/holder)*np.min(pvals))

def get_args():
	#get the arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('tar',type=int)
	parser.add_argument('prsize',type=int)
	parser.add_argument('pval',type=float)
	args = parser.parse_args()
	return args

def main():
	#let the command line arguments flow in
	args = get_args()
	
	#import data. same old pandas into non pandas trick, now with a "header"
	#(the first line of the file is pretty useless)
	fimo = pd.read_csv('fimo.txt', sep='\t', index_col=None, header=0)
	fimo = fimo.values
	#just check if the file is empty. some motifs hit nothing
	if not fimo.size:
		sys.exit(0)

	#let's kick out the potential overlapping motif instances
	del_inds = []
	prev_ind = 0
	motifs = []
	for i in range(1,fimo.shape[0]):
		#case number one - the motifs are for the same gene
		if (fimo[i,0] == fimo[prev_ind,0] and fimo[i,1] == fimo[prev_ind,1]):
			#so we need to check if we gots ourselves an overlap
			if overlapcheck(fimo[i,[2,3]],fimo[prev_ind,[2,3]]):
				#if we're in here, we have an overlap
				if fimo[i,6] < fimo[prev_ind,6]:
					#if we're in here, our current motif outranks the previous one on p-value
					del_inds.append(prev_ind)
					prev_ind = i
				else:
					#if we're in here, the previous one's p-value was stronger
					del_inds.append(i)
			else:
				#if we're in here, then there's no overlap
				prev_ind = i
		else:
			#case number two - different motif instance altogether
			#is it a new motif?
			if fimo[i,0] != fimo[prev_ind,0]:
				motifs.append(fimo[prev_ind,0])
			#shift index
			prev_ind = i
	#we're not actually including the final unique motif in motifs
	motifs.append(fimo[prev_ind,0])
	#so, the digging for redundant motifs has completed. let's see what we found
	fimo = np.delete(fimo,del_inds,0)
	
	#how many tests will we actually have?
	test_count = []
	for mot in motifs:
		test_count.append(len(np.unique(fimo[fimo[:,0]==mot,1])))
	test_count = np.asarray(test_count)
	
	#time to get testing
	writer = open('fimo_found.txt','a')
	start_pos = 0
	test_counter = 0
	for i in range(1,fimo.shape[0]):
		if (fimo[i,0] != fimo[start_pos,0] or fimo[i,1] != fimo[start_pos,1]):
			#we've got a motif flip going. test from start_pos to i-1
			pvalues = fimo[start_pos:i,6]
			pvalues = np.sort(pvalues)
			binom_p = []
			#how many "flips" will we be making
			motsize = fimo[start_pos,3] - fimo[start_pos,2] + 1
			flips = 2 * (args.prsize - motsize + 1)
			#test the top N motifs
			#obviously if we have less than N, test what we do have
			for j in range(np.min([(i-start_pos),args.tar])):
				#this is our binomial probability
				coin = hommel(pvalues[:(j+1)])
				#the actual binomial test is done using .sf
				#but you have to go one less than the actual number of successes you expect
				#as .sf is 1-CDF, so including the number of successes you wants gets it subtracted away
				binom_p.append(sps.binom.sf(j,flips,coin))
			#convert to array for ease of things
			binom_p = np.asarray(binom_p)
			#perform FDR (good old Bonferroni)
			binom_p = binom_p * test_count[test_counter]
			#do we have a hit?
			if any(binom_p <= args.pval):
				writer.write(fimo[start_pos,0]+'\t'+fimo[start_pos,1]+'\n')
			#well, that's a wrap. set it up for the next one
			#potentially wedge forward the Bonferroni counter
			if fimo[i,0] != fimo[start_pos,0]:
				test_counter = test_counter+1
			#definitely wedge forward the hit start
			start_pos = i
	#however, at the end, we need to process the last one
	pvalues = fimo[start_pos:,6]
	pvalues = np.sort(pvalues)
	binom_p = []
	#how many "flips" will we be making
	motsize = fimo[start_pos,3] - fimo[start_pos,2] + 1
	flips = 2 * (args.prsize - motsize + 1)
	#test the top N motifs
	#obviously if we have less than N, test what we do have
	for j in range(np.min([(i-start_pos),args.tar])):
		#this is our binomial probability
		coin = hommel(pvalues[:(j+1)])
		#the actual binomial test is done using .sf
		#but you have to go one less than the actual number of successes you expect
		#as .sf is 1-CDF, so including the number of successes you wants gets it subtracted away
		binom_p.append(sps.binom.sf(j,flips,coin))
	#convert to array for ease of things
	binom_p = np.asarray(binom_p)
	#perform FDR (good old Bonferroni)
	binom_p = binom_p * test_count[test_counter]
	#do we have a hit?
	if any(binom_p <= args.pval):
		writer.write(fimo[start_pos,0]+'\t'+fimo[start_pos,1]+'\n')
	#close the writer
	writer.close()

#Sam trick again
if __name__ == "__main__":
	main()