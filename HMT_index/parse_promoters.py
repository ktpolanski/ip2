import numpy as np
import pandas as pd

reader = open('promoters_rough.fa','r')
bedfile = pd.read_csv('promoters.bed', sep='\t', index_col=None, header=None)
bedfile = bedfile.values
bedfile = bedfile[:,3]
bedfile = np.asarray(['>'+bed+'\n' for bed in bedfile])
fafile = np.asarray(reader.readlines())
reader.close()
fafile[np.arange(0,len(fafile),2)] = bedfile
writer = open('promoters.fa','w')
writer.writelines(fafile)
writer.close()