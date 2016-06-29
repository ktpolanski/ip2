import pandas as pd
import numpy as np

prom = pd.read_csv('promoters.bed',header=None,index_col=None,sep='\t').values
fid = open('promoter_lengths.txt','w')
toggle = ''

for i in np.arange(prom.shape[0]):
	fid.write(toggle+prom[i,3]+'\t'+str(prom[i,2]-prom[i,1]))
	toggle='\n'

fid.close()