import numpy as np
import pandas as pd
import argparse

def main():
	#sponge up the input
	parser = argparse.ArgumentParser()
	parser.add_argument('--Input', dest='input', type=argparse.FileType('r'), nargs='+', required=True, help='The CSV marginal files produced by hCSI (or oCSI)')
	parser.add_argument('--Threshold', dest='thresh', type=float, required=True, help='The desired probability threshold for making the network binary')
	args = parser.parse_args()
	#process each file independently
	for i in np.arange(len(args.input)):
		inp = pd.read_csv(args.input[i],dtype=str,index_col=0,header=0)
		inp = inp.astype(float)
		#zero the edges that are not strong enough
		inp[inp<args.thresh] = 0
		#kick out the irrelevant genes, i.e. those that are not parents or children to anything
		#(if they aren't, then the sum of their marginal row and column will be 0)
		sums = np.sum(inp,axis=0) + np.sum(inp,axis=1)
		keep = np.where(sums>0)[0]
		inp = inp.iloc[keep,keep]
		#now the marginal matrix is sufficiently filtered and we can get exporting
		writer = open(args.input[i].name+'.gml','w')
		#if we have an empty matrix as a result of too stringent a threshold, just make an empty file
		if not inp.empty:
			#header
			writer.write('graph [\ndirected 1\n')
			#list of nodes
			nodes = list(inp.index)
			for j in np.arange(len(nodes)):
				writer.write('node [\nid '+str(j)+'\nlabel "'+nodes[j]+'"\n]\n')
			#list of edges
			for j in np.arange(len(nodes)):
				for k in np.arange(len(nodes)):
					if inp.iloc[j,k]>0:
						#we have an edge
						writer.write('edge [\nsource '+str(k)+'\ntarget '+str(j)+'\nweight '+str(inp.iloc[j,k])+'\ninteraction "pd"\n]\n')
			#close the thing off
			writer.write(']')
		writer.close()

if __name__ == '__main__':
    main()