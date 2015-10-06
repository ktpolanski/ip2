import argparse
import sys
import numpy as np
import pandas as pd
import wigwams as ww
import time

def parse_dfs():
	'''
	Take all things command line and data and return nice, handy, PANDAS parsed versions of gene expression, DEGs and other arguments.
	'''
	
	#accept all useful arguments in the command line, most of them optional
	#literally the only required one being the expression CSV, as that's just needed for Wigwams to mine
	#the rest has defaults in place, as used during the PRESTA runs
	parser = argparse.ArgumentParser()
	parser.add_argument('--Expression', dest='expr', type=argparse.FileType('r'), required=True, help='CSV file featuring the expression profiles of the genes across all the conditions. Gene names (case insensitive) in first column. First two rows for headers - row one for condition name, row two for time point')
	parser.add_argument('--DEGs', dest='deg', type=argparse.FileType('r'), help='Optional binary CSV file featuring differential expression information across all the conditions. Gene names (case insensitive) in first column. Condition names in first row.')
	parser.add_argument('--PoolNumber', dest='pool', type=int, default=1, help='Number of processes to run in module mining (to potentially parallelise the mining and accelerate the code). Default: 1')
	parser.add_argument('--Mining_SetSizes', dest='sets', type=int, default=[50, 100, 150, 200, 250], nargs='+', help='Module mining. Top co-expressed gene set sizes to be scanned for evidence of dependent co-expression. Provide as space delimited list. Default: 50 100 150 200 250')
	parser.add_argument('--Mining_Alpha', dest='alpha', default=0.05, type=float, help='Module mining. The significance threshold for dependent co-expression testing, is Bonferroni corrected within the script. Default: 0.05')
	parser.add_argument('--Mining_CorrelationNet', dest='corrnet', default=0.7, type=float, help='Module mining. In the case of scarce profiles, dependent co-expression testing is stopped if the top correlated genes stop being at least this PCC-correlated with the seed. Default: 0.7')
	parser.add_argument('--Merging_Overlap', dest='merging_overlap', default=0.3, type=float, help='Module merging. The fraction of a module\'s members that have to be present in the overlap with another module spanning the same conditions to trigger merging. Default: 0.3')
	parser.add_argument('--Merging_MeanCorrelation', dest='meancorr', default=0.9, type=float, help='Module merging. Initialise merging if the mean expression profiles of two modules across all of the conditions they span are at least this PCC-correlated, regardless of membership overlap. Default: 0.9')
	parser.add_argument('--Merging_CorrelationFilter', dest='corrfilt', default=0.8, type=float, help='Module merging. If a gene\'s profile isn\'t at least this PCC-correlated with the mean profile of the larger module, don\'t transfer it over when merging, discarding it instead. Default: 0.8')
	parser.add_argument('--Sweeping_Overlap', dest='sweep_overlap', default=0.5, type=float, help='Module sweeping. Discard a module spanning a subset of another module\'s conditions if its membership is made up of at least this proportion of the module with the condition superset. Default: 0.5')
	parser.add_argument('--SizeThresholds', dest='thresh', default=[], type=int, nargs='+', help='Module size thresholding. Optional post-processing step to filter down the module list to modules of at least a given size. Provide as space delimited list, with first element being the desired minimal size of 2-condition modules and the last being the desired minimal size of modules across all the conditions. Default: [] (procedure off)')
	args = parser.parse_args()

	#parse the expression CSV
	print('Parsing expression CSV file...')
	expr_df = pd.read_csv(args.expr, index_col=0, header=[0,1])
	#case insensitive gene names = caps all the gene names on sight
	expr_df.index = map(lambda x:x.upper(), expr_df.index.tolist())
	#currently the time points are stored as strings, which is not perfect
	#make a new bunch of tuples
	old_colnames = expr_df.columns.values.tolist()
	new_colnames = []
	for lab in old_colnames:
		new_colnames.append((lab[0].title(), float(lab[1])))
	#the names=[] part allows for easy unearthing of all things treatment for correlations later
	expr_df.columns = pd.MultiIndex.from_tuples(new_colnames,names=['Treatment','Time'])

	#parse the DEG CSV
	if args.deg is not None:
		print('Parsing DEG CSV file...')
		deg_df = pd.read_csv(args.deg, index_col=0, header=0)
		#once again, caps gene names on sight
		deg_df.index = map(lambda x:x.upper(), deg_df.index.tolist())
		deg_df.columns = map(lambda x:x.title(), deg_df.columns.tolist())
		
		#filter to the same genes in both data frames
		genes = set(deg_df.index) & set(expr_df.index)
		deg_df = deg_df[deg_df.index.isin(genes)]
		expr_df = expr_df[expr_df.index.isin(genes)]
		
		#filter to the same columns in both data frames
		#need to walk around .loc[] not working with this for reasons I don't understand
		expconds = expr_df.columns.get_level_values('Treatment').tolist()
		degconds = deg_df.columns.tolist()
		intconds = list(set(expconds) & set(degconds))
		#filtering involves using .iloc[] as that somehow works
		#dig up the indices of the condition overlap
		expmask = []
		for i in range(0,len(expconds)):
			if expconds[i] in intconds:
				expmask.append(i)
		degmask = []
		for i in range(0,len(degconds)):
			if degconds[i] in intconds:
				degmask.append(i)
		#and now that we have the two sets of indices, filter the columns
		expr_df = expr_df.iloc[:,expmask]
		deg_df = deg_df.iloc[:,degmask]
		
		#now, filter out genes that are insufficiently DE
		deg_df = deg_df[deg_df.sum(axis=1)>0]
		genes = deg_df.index.tolist()
		expr_df = expr_df[expr_df.index.isin(genes)]
	else:
		#just make a mock one if there is none, makes it easier for the script later on
		#no need for filtering as obviously all is well
		print('Creating mock DEG dataframe...')
		expgenes = expr_df.index.tolist()
		expconds = list(set(expr_df.columns.get_level_values('Treatment')))
		deg_df = pd.DataFrame(np.ones((len(expgenes),len(expconds))), index=expgenes, columns=expconds)

	#time for a few sanity checks. need at least 2 conditions (and preferably more than 2 genes too)
	if 0 in deg_df.shape or 1 in deg_df.shape:
		print('error: insufficient number of conditions and/or genes')
		sys.exit()
	
	#need some genes to mine as well
	if all(np.sum(deg_df.values,axis=1)<2):
		print('error: no genes differentially expressed across two or more conditions')
		sys.exit()
	
	#also would be nice that if we have the thresholding sizes, they actually match in number to condition span
	if args.thresh:
		if len(args.thresh) != (len(deg_df.columns)-1):
			print('error: incompatible number of thresholding set sizes and actual number of data sets')
			sys.exit()

	#time for data standardisation
	#temporarily commented out to evaluate run
	'''
	print('Standardising data...')
	for cond in list(set(expr_df.columns.get_level_values('Treatment'))):
		vals = expr_df.loc[:,cond].values
		expr_df.loc[:,cond] = (vals - np.mean(vals,axis=1)[:,None]) / np.std(vals,axis=1)[:,None]
	'''
	
	#I think we have everything this function can do. spit it out for further consumption
	return (expr_df, deg_df, args)

def wigwams_analysis_default(expr_df, deg_df, args):
	'''
	The bog standard Wigwams pipeline, as described in the paper: mining -> merging -> sweeping
	
	Input:
		* expr_df - PANDAS data frame of gene expression, as parsed by parse_dfs(). Multiindexing of columns - treatments and times.
		* deg_df - PANDAS data frame of DEG status, as parsed by parse_dfs().
		* args - command line arguments, as parsed by parse_dfs().
	'''
	#mining. the time consuming part
	t0 = time.time()
	ww.mining(expr_df, deg_df, pool=args.pool, sets=args.sets, alpha = args.alpha, corrnet=args.corrnet)
	t1 = time.time()
	m, s = divmod(round(t1-t0), 60)
	h, m = divmod(m, 60)
	print('Took %d:%02d:%02d. Module mining complete.' % (h, m, s))
	#merging
	ww.merging(expr_df, overlap=args.merging_overlap, meancorr=args.meancorr, corrfilt=args.corrfilt)

def main():
	#parse the inputs
	(expr_df, deg_df, args) = parse_dfs()
	#default analysis pipeline. split off into function to experiment with workflows if desired
	wigwams_analysis_default(expr_df, deg_df, args)

#this is a pretty cool Sam trick
#if I put the whole code in main() then I can access it from something like iPython
#and at the same time have more functions in here if I so desire
#whilst this little bit here at the bottom makes it so that command line fires main()
if __name__ == "__main__":
	main()