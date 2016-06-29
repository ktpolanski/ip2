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

def bh(pvalues):
	'''
	Computes the Benjamini-Hochberg FDR correction.
	
	Input:
		* pvals - vector of p-values to correct
	'''
	pvalues = np.array(pvalues)
	n = float(pvalues.shape[0])
	new_pvalues = np.empty(n)
	values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
	values.sort()
	values.reverse()
	new_values = []
	for i, vals in enumerate(values):
		rank = n - i
		pvalue, index = vals
		new_values.append((n/rank) * pvalue)
	for i in range(0, int(n)-1):
		if new_values[i] < new_values[i+1]:
			new_values[i+1] = new_values[i]
	for i, vals in enumerate(values):
		pvalue, index = vals
		new_pvalues[index] = new_values[i]
	return new_pvalues

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--FoundMotifs', dest='motifs', type=argparse.FileType('r'), required=True, help='HMT_index output fimo_found.txt file, which features instances of motifs identified as present in genes\' promoters.')
	parser.add_argument('--Universe', dest='universe', type=argparse.FileType('r'), required=True, help='TXT file with the list of all gene IDs present at the HMT_index step, one gene ID per line.')
	parser.add_argument('--Input', dest='input', type=argparse.FileType('r'), required=True, help='TXT file with MEME-LaB-friendly parsing of gene groups to analyse (tab delimited, column one group ID, column two gene ID.')
	parser.add_argument('--MotifAnnotation', dest='mot_annot', default=None, type=argparse.FileType('r'), help='CSV file with additional information on the motifs to include in the export.')
	parser.add_argument('--Alpha', dest='alpha', default=0.05, type=float, help='Desired significance threshold. FDR correction is performed both with the Bonferroni and Benjamini-Hochberg methods within the script. Default: 0.05')
	args = parser.parse_args()
	return args

def main():
	args = parse_args()
	#import stuff
	foundmotifs = pd.read_csv(args.motifs, sep='\t', index_col=None, header=None)
	foundmotifs = foundmotifs.values
	#upper the foundmotifs as well
	for i in np.arange(foundmotifs.shape[0]):
		for j in np.arange(foundmotifs.shape[1]):
			foundmotifs[i,j] = foundmotifs[i,j].upper()
	motifs = foundmotifs[:,0]
	motifs = np.array([item.upper() for item in motifs])
	motifs = np.unique(motifs)
	#no need to open up a reader for the universe, argparse does it automatically
	universe = args.universe.read().splitlines()
	universe = np.array([item.upper() for item in universe])
	universe = np.unique(universe)
	input = pd.read_csv(args.input, sep='\t', index_col=None, header=None)
	input = input.values
	input[:,1] = np.array([item.upper() for item in input[:,1]])
	groups = np.unique(input[:,0])
	strgroups = [str(item) for item in groups]
	#optional annotation input, and check if it works
	mot_annot = None
	if args.mot_annot:
		mot_annot = pd.read_csv(args.mot_annot, sep=',', index_col=None, header=None)
		mot_annot = mot_annot.values
		mot_annot[:,0] = np.array([item.upper() for item in mot_annot[:,0]])
		if len(np.unique(mot_annot[:,0])) != mot_annot.shape[0]:
			sys.stdout.write('Duplicate motifs present in annotation. Ignoring annotation')
			mot_annot = None
		else:
			for mot in motifs:
				if mot not in mot_annot[:,0]:
					sys.stdout.write('Failed to identify motif '+mot+' in the provided annotation. Ignoring annotation.\n')
					mot_annot = None
					break
	#if there's still a mot_annot, it includes the motif space of motifs (the variable)
	motifs.sort()
	if mot_annot is not None:
		#safekeeping the header
		header = mot_annot[0,:]
		mot_annot = np.delete(mot_annot,0,0)
		#wiping the non-motifs motifs
		del_inds = []
		for i in range(mot_annot.shape[0]):
			if mot_annot[i,0] not in motifs:
				del_inds.append(i)
		if del_inds:
			mot_annot = np.delete(mot_annot,del_inds,0)
		#sorting to be the same order
		mask = np.asarray([i[0] for i in sorted(enumerate(mot_annot[:,0]), key=lambda x:x[1])])
		mot_annot = mot_annot[mask,:]
	#note the lack of an else - it's entirely possible that the mot_annot gets None'd above
	if mot_annot is None:
		#make a "mock annot" to have the export be homogeneous
		header = ['Motif ID']
		mot_annot = motifs[:,None]
	#pre-generate p-value storage space
	pvalues = np.zeros((len(motifs),len(groups)))
	pvalues_bonf = np.zeros((len(motifs),len(groups)))
	pvalues_bh = np.zeros((len(motifs),len(groups)))
	pvalues_glbonf = np.zeros((len(motifs),len(groups)))
	#pre-generate overrepresentation holders
	over_head = list(header)
	over_head.extend(['Module ID', 'P-Value', 'Bound Genes'])
	over_bonf = open('Overrepresentation_Bonferroni.txt','w')
	over_bh = open('Overrepresentation_Benjamini-Hochberg.txt','w')
	over_glbonf = open('Overrepresentation_GlobalBonferroni.txt','w')
	over_bonf.write('\t'.join(over_head))
	over_bh.write('\t'.join(over_head))
	over_glbonf.write('\t'.join(over_head))
	#dictionaries are efficient
	matdict = {}
	for mot in mot_annot[:,0]:
		motgenes = np.unique(foundmotifs[np.where(foundmotifs[:,0]==mot)[0],1])
		motgenes = np.array(list(set(motgenes).intersection(set(universe))))
		matdict[mot] = motgenes
	#testing proper. loop over gene groups
	for i in range(len(groups)):
		genes = input[np.where(input[:,0]==groups[i])[0],1]
		genes = np.array(list(set(genes).intersection(set(universe))))
		for j in range(len(motifs)):
			motgenes = matdict[motifs[j]]
			hitgenes = np.array(list(set(genes).intersection(set(motgenes))))
			#differences with Matlab version are here
			#due to the fact we're using a wider universe
			pvalues[j,i] = np.exp(log_hypergeometric(len(hitgenes),len(motgenes),len(genes),len(universe)))
			if pvalues[j,i]>1:
				pvalues[j,i]=1
		#"local" bonferroni - just correct within a group
		pvalues_bonf[:,i] = pvalues[:,i] * len(motifs)
		#benjamini-hochberg is also within a group
		pvalues_bh[:,i] = bh(pvalues[:,i])
		#"global" bonferroni corrects across groups as well
		pvalues_glbonf[:,i] = pvalues[:,i] * len(groups) * len(motifs)
		for j in range(len(motifs)):
			if pvalues_bonf[j,i]>1:
				pvalues_bonf[j,i]=1 
			elif pvalues_bonf[j,i] <= args.alpha:
				motgenes = matdict[motifs[j]]
				hitgenes = np.array(list(set(genes).intersection(set(motgenes))))
				hitgenes.sort()
				line = list(mot_annot[j,:])
				line.extend([strgroups[i], str(pvalues_bonf[j,i]), ';'.join(hitgenes)])
				over_bonf.write('\n'+'\t'.join(line))
			if pvalues_bh[j,i]>1:
				pvalues_bh[j,i]=1
			elif pvalues_bh[j,i] <= args.alpha:
				motgenes = matdict[motifs[j]]
				hitgenes = np.array(list(set(genes).intersection(set(motgenes))))
				hitgenes.sort()
				line = list(mot_annot[j,:])
				line.extend([strgroups[i], str(pvalues_bh[j,i]), ';'.join(hitgenes)])
				over_bh.write('\n'+'\t'.join(line))
			if pvalues_glbonf[j,i]>1:
				pvalues_glbonf[j,i]=1 
			elif pvalues_glbonf[j,i] <= args.alpha:
				motgenes = matdict[motifs[j]]
				hitgenes = np.array(list(set(genes).intersection(set(motgenes))))
				hitgenes.sort()
				line = list(mot_annot[j,:])
				line.extend([strgroups[i], str(pvalues_glbonf[j,i]), ';'.join(hitgenes)])
				over_glbonf.write('\n'+'\t'.join(line))
	#time to finish off the exports
	writer1 = open('Full_Raw_P-Values.txt','w')
	writer2 = open('Full_Bonferroni_P-Values.txt','w')
	writer3 = open('Full_Benjamini-Hochberg_P-Values.txt','w')
	writer4 = open('Full_GlobalBonferroni_P-Values.txt','w')
	line = list(header)
	line.extend(strgroups)
	writer1.write('\t'.join(line))
	writer2.write('\t'.join(line))
	writer3.write('\t'.join(line))
	writer4.write('\t'.join(line))
	#now, the actual lines
	for i in range(pvalues.shape[0]):
		line1 = list(mot_annot[i,:])
		line1.extend([str(item) for item in pvalues[i,:]])
		line2 = list(mot_annot[i,:])
		line2.extend([str(item) for item in pvalues_bonf[i,:]])
		line3 = list(mot_annot[i,:])
		line3.extend([str(item) for item in pvalues_bh[i,:]])
		line4 = list(mot_annot[i,:])
		line4.extend([str(item) for item in pvalues_glbonf[i,:]])
		writer1.write('\n'+'\t'.join(line1))
		writer2.write('\n'+'\t'.join(line2))
		writer3.write('\n'+'\t'.join(line3))
		writer4.write('\n'+'\t'.join(line4))
	writer1.close()
	writer2.close()
	writer3.close()
	writer4.close()
	#now, time for the significant aspect
	mask_bonf1 = np.nonzero(np.sum(pvalues_bonf<=args.alpha,0))[0]
	mask_bonf2 = np.nonzero(np.sum(pvalues_bonf<=args.alpha,1))[0]
	mask_bh1 = np.nonzero(np.sum(pvalues_bh<=args.alpha,0))[0]
	mask_bh2 = np.nonzero(np.sum(pvalues_bh<=args.alpha,1))[0]
	mask_glbonf1 = np.nonzero(np.sum(pvalues_glbonf<=args.alpha,0))[0]
	mask_glbonf2 = np.nonzero(np.sum(pvalues_glbonf<=args.alpha,1))[0]
	strgroups_bonf = np.array(strgroups)[mask_bonf1]
	strgroups_bh = np.array(strgroups)[mask_bh1]
	strgroups_glbonf = np.array(strgroups)[mask_glbonf1]
	mot_annot_bonf = mot_annot[mask_bonf2,:]
	mot_annot_bh = mot_annot[mask_bh2,:]
	mot_annot_glbonf = mot_annot[mask_glbonf2,:]
	pvals_bonf = pvalues_bonf[mask_bonf2[:, None], mask_bonf1]
	pvals_bh = pvalues_bh[mask_bh2[:, None], mask_bh1]
	pvals_glbonf = pvalues_glbonf[mask_glbonf2[:, None], mask_glbonf1]
	writer1 = open('Significant_Bonferroni_P-Values.txt','w')
	writer2 = open('Significant_Benjamini-Hochberg_P-Values.txt','w')
	writer3 = open('Significant_GlobalBonferroni_P-Values.txt','w')
	line = list(header)
	line.extend(list(strgroups_bonf))
	writer1.write('\t'.join(line))
	line = list(header)
	line.extend(list(strgroups_bh))
	writer2.write('\t'.join(line))
	line = list(header)
	line.extend(list(strgroups_glbonf))
	writer3.write('\t'.join(line))
	for i in range(pvals_bonf.shape[0]):
		line = list(mot_annot_bonf[i,:])
		line.extend([str(item) for item in pvals_bonf[i,:]])
		writer1.write('\n'+'\t'.join(line))
	for i in range(pvals_bh.shape[0]):
		line = list(mot_annot_bh[i,:])
		line.extend([str(item) for item in pvals_bh[i,:]])
		writer2.write('\n'+'\t'.join(line))
	for i in range(pvals_glbonf.shape[0]):
		line = list(mot_annot_glbonf[i,:])
		line.extend([str(item) for item in pvals_glbonf[i,:]])
		writer3.write('\n'+'\t'.join(line))
	writer1.close()
	writer2.close()
	writer3.close()
	#generate output for that webapp thing
	mot_annot = motifs[:,None]
	writer1 = open('Webapp_p_values.txt','w')
	writer2 = open('Webapp_table_gene_counts.txt','w');
	writer3 = open('Webapp_which_genes.txt','w');
	writer4 = open('Webapp_clusters.txt','w');
	writer1.write('\t'+'\t'.join(strgroups)+'\n')
	writer2.write('Motif\tUniverse\t'+'\t'.join(strgroups)+'\n\t'+str(len(universe)))
	writer3.write('Motif\tFrequency in Universe\tCluster\tCluster Size\tFrequency in Cluster\tGene IDs\tp-value for over-representation\n')
	for i in range(len(groups)):
		genes = input[np.where(input[:,0]==groups[i])[0],1]
		genes = np.array(list(set(genes).intersection(set(universe))))
		writer1.write('\t'+str(len(genes)))
		writer2.write('\t'+str(len(genes)))
		for j in range(len(genes)):
			writer4.write(strgroups[i]+'\t'+genes[j]+'\n')
	for j in range(pvalues.shape[0]):
		line = list(mot_annot[j,:])
		line.extend([str(item) for item in pvalues[j,:]])
		writer1.write('\n'+'\t'.join(line))
		motgenes = matdict[motifs[j]]
		writer2.write('\n'+motifs[j]+'\t'+str(len(motgenes)))
		for i in range(pvalues.shape[1]):
			genes = input[np.where(input[:,0]==groups[i])[0],1]
			genes = np.array(list(set(genes).intersection(set(universe))))
			hitgenes = np.array(list(set(genes).intersection(set(motgenes))))
			hitgenes.sort()
			writer2.write('\t'+str(len(hitgenes)))
			writer3.write('\n'+motifs[j]+'\t'+str(len(motgenes))+'\t'+strgroups[i]+'\t'+str(len(genes))+'\t'+str(len(hitgenes))+'\t'+','.join(hitgenes)+'\t'+str(pvalues[j,i]))
	writer1.close()
	writer2.close()
	writer3.close()
	writer4.close()

if __name__ == "__main__":
	main()