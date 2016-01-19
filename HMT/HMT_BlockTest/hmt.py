import argparse
import numpy as np
import pandas as pd
import sys

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

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--FoundMotifs', dest='motifs', type=argparse.FileType('r'), required=True, help='HMT_index output fimo_found.txt file, which features instances of motifs identified as present in genes\' promoters.')
	parser.add_argument('--Universe', dest='universe', type=argparse.FileType('r'), required=True, help='TXT file with the list of all gene IDs present at the HMT_index step, one gene ID per line.')
	parser.add_argument('--BlockSize', dest='block', default=7, type=int, help='How many adjacent neighbours to evaluate.')
	parser.add_argument('--Alpha', dest='alpha', default=0.05, type=float, help='Desired significance threshold. FDR correction is performed both with the Bonferroni and Benjamini-Hochberg methods within the script. Default: 0.05')
	args = parser.parse_args()
	return args

def main():
	args = parse_args()
	#import stuff
	foundmotifs = pd.read_csv(args.motifs, sep='\t', index_col=None, header=None)
	foundmotifs = foundmotifs.values
	motifs = foundmotifs[:,0]
	motifs = np.array([item.upper() for item in motifs])
	motifs = np.unique(motifs)
	#no need to open up a reader for the universe, argparse does it automatically
	universe = args.universe.read().splitlines()
	universe = np.array([item.upper() for item in universe])
	#if there's still a mot_annot, it includes the motif space of motifs (the variable)
	motifs.sort()
	#make a "mock annot" to have the export be homogeneous
	header = ['Motif ID']
	mot_annot = motifs[:,None]
	#pre-generate p-value storage space
	pvalues = np.zeros((len(motifs),1))
	pvalues_bonf = np.zeros((len(motifs),1))
	#dictionaries are efficient
	matdict = {}
	for mot in mot_annot[:,0]:
		motgenes = np.unique(foundmotifs[np.where(foundmotifs[:,0]==mot)[0],1])
		motgenes = np.array(list(set(motgenes).intersection(set(universe))))
		matdict[mot] = motgenes
	#peanut butter export time
	writer = open('Found.txt','w')
	#testing proper. loop over gene groups
	for i in range(len(universe)-args.block+1):
		genes = universe[i:(i+args.block)]
		#check if they're all on the same chromosome - first four letters will be the same
		gold = genes[0][0:4]
		contflip = False
		for j in range(1,len(genes)):
			#we have something on a different chromosome, ladies and germs
			if not gold==genes[j][0:4]:
				contflip = True
		if contflip:
			continue
		for j in range(len(motifs)):
			motgenes = matdict[motifs[j]]
			hitgenes = np.array(list(set(genes).intersection(set(motgenes))))
			#differences with Matlab version are here
			#due to the fact we're using a wider universe
			pvalues[j] = np.exp(log_hypergeometric(len(hitgenes),len(motgenes),len(genes),len(universe)))
			if pvalues[j]>1:
				pvalues[j]=1
		pvalues_bonf = pvalues * len(motifs)
		pvalues_bonf[pvalues_bonf>1] = 1
		writer.write(universe[i]+'\t'+str(np.min(pvalues_bonf))+'\n')
		
if __name__ == "__main__":
	main()