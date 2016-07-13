import pandas as pd
import numpy as np

#soak up current promoters.bed file
prom = pd.read_csv('promoters.bed',header=None,index_col=None,sep='\t').values
annot = pd.read_csv('annot.gff3',header=None,index_col=None,sep='\t',comment='#').values

#the order is the same in the promoters.bed as it is in annot.gff3
#after all, annot.gff3 was used to craft promoters.bed
geneinds = []
for i in np.arange(annot.shape[0]):
	if 'gene' in annot[i,:]:
		geneinds.append(i)
#and need an "end of file" geneinds as well
geneinds.append(annot.shape[0]+1)

#loop over the genes
for i in np.arange(prom.shape[0]):
	#dig out the subannot and filter it to cds
	subannot = annot[geneinds[i]:geneinds[i+1],:]
	noncds = []
	for j in np.arange(subannot.shape[0]):
		if 'CDS' not in subannot[j,:]:
			noncds.append(j)
	temp = np.delete(subannot,noncds,axis=0)
	#if there's anything found, continue
	if temp.size:
		if prom[i,5]=='+':
			#"left to right" transcript. find earliest CDS start and replace end of promoter
			newpos = np.min(temp[:,3])
			prom[i,2] = newpos
		else:
			#"right to left" transcript. find latest CDS "end", which is in fact the earliest CDS start
			#and replace "start" of promoter which is in fact its end
			newpos = np.max(temp[:,4])
			prom[i,1] = newpos

#export the 5' UTR'd promoter file
np.savetxt('promoters.bed', prom, fmt='%s\t%i\t%i\t%s\t%i\t%s')