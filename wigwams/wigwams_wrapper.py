import matplotlib
matplotlib.use('Agg')
import argparse
import numpy as np
import pandas as pd
import wigwams as ww
import os
import sys

def parse_args():
	'''
	Parse all things command line, and then return them
	'''
	
	#accept all useful arguments in the command line, most of them optional
	#literally the only required one being the expression CSV, as that's just needed for Wigwams to mine
	#the rest has defaults in place, as used during the PRESTA runs
	parser = argparse.ArgumentParser()
	parser.add_argument('--Expression', dest='expr', type=argparse.FileType('r'), required=True, help='CSV file featuring the expression profiles of the genes across all the conditions. Gene names (case insensitive) in first column. First two rows for headers - row one for condition name, row two for time point.')
	parser.add_argument('--DEGs', dest='deg', default=None, type=argparse.FileType('r'), help='Optional binary CSV file featuring differential expression information across all the conditions. Gene names (case insensitive) in first column. Condition names in first row.')
	parser.add_argument('--NoStandardising', dest='stand', action='store_false', help='Flag. If provided, the data isn\'t standardised (scaled to N(0,1)) on input.')
	parser.add_argument('--Legacy', dest='leg', action='store_true', help='Flag. If provided, runs algorithm in legacy mode (as described in Polanski et al 2014) in terms of module mining and merging. If not provided, slightly improves those procedures.')
	parser.add_argument('--NonRedundantOutput', dest='nro', action='store_true', help='Flag. If provided, runs algorithm in iterational redundancy removal mode which yields considerably smaller and less redundant output, but at the cost of some information loss.')
	parser.add_argument('--PoolNumber', dest='pool', type=int, default=1, help='Number of processes to run in module mining (to potentially parallelise the mining and accelerate the code). Default: 1 (not parallel)')
	parser.add_argument('--Mining_SetSizes', dest='sets', type=int, default=[50, 100, 150, 200, 250], nargs='+', help='Module mining. Top co-expressed gene set sizes to be scanned for evidence of dependent co-expression. Provide as space delimited list. Default: 50 100 150 200 250')
	parser.add_argument('--Mining_Alpha', dest='alpha', default=0.05, type=float, help='Module mining. The significance threshold for dependent co-expression testing, is Bonferroni corrected within the script. Default: 0.05')
	parser.add_argument('--Mining_CorrelationNet', dest='corrnet', default=0.7, type=float, help='Module mining. In the case of scarce profiles, dependent co-expression testing is stopped if the top correlated genes stop being at least this PCC-correlated with the seed. Default: 0.7')
	parser.add_argument('--Merging_Overlap', dest='merging_overlap', default=0.3, type=float, help='Module merging. The fraction of a module\'s members that have to be present in the overlap with another module spanning the same conditions to trigger merging. Default: 0.3')
	parser.add_argument('--Merging_MeanCorrelation', dest='meancorr', default=0.9, type=float, help='Module merging. Initialise merging if the mean expression profiles of two modules across all of the conditions they span are at least this PCC-correlated, regardless of membership overlap. Default: 0.9')
	parser.add_argument('--Merging_CorrelationFilter', dest='corrfilt', default=0.8, type=float, help='Module merging. If a gene\'s profile isn\'t at least this PCC-correlated with the mean profile of the larger module, don\'t transfer it over when merging, discarding it instead. Default: 0.8')
	parser.add_argument('--Sweeping_Overlap', dest='sweeping_overlap', default=0.5, type=float, help='Module sweeping. Discard a module spanning a subset of another module\'s conditions if its membership is made up of at least this proportion of the module with the condition superset. In non redundant output mode this becomes the upper bound of the sweeping thresholds to evaluate. Default: 0.5')
	parser.add_argument('--Sweeping_Overlap_LowBound', dest = 'sweep_low', default = 0.1, type=float, help='Module sweeping. In non redundant output mode, the lower bound of sweeping to evaluate, with 0.05 resolution between this and the upper bound. Ignored in legacy mode. Default: 0.1')
	parser.add_argument('--SizeThresholds', dest='thresh', default=None, type=int, nargs='+', help='Module size thresholding. Optional post-processing step to filter down the module list to modules of at least a given size. Provide as space delimited list, with first element being the desired minimal size of 2-condition modules and the last being the desired minimal size of modules across all the conditions. Default: None (procedure off)')
	parser.add_argument('--Export_Annotation', dest='annot', default=None, type=argparse.FileType('r'), help='TSV (tab-separated file) with an annotation carrying additional information on the genes in the mining to include in the final export. First column must match the identifiers used, second column must be a public identifier, third column onwards optional. Please remove any "unmapped" information from the mapping. Default: None (no additional annotating of export)')
	parser.add_argument('--Export_Hyperlink', dest='hyper', default=None, type=str, help='An optional extension to the annotation, adding Excel-friendly hyperlinks to genes identified in modules to the module export in a particular online resource. Include a generic link to finding the gene IDs in the second column of the annotation in the online resource of choice, with the second column identifier place indicated with {gene}. Please provide wrapped in quotes. Default: None (no additional hyperlinks)')
	parser.add_argument('--JobName', dest='job', default='job', type=str, help='Job name for output naming purposes. Default: job')
	args = parser.parse_args()
	return args

def parse_expr_df(args):
	'''
	Parse the expression CSV within parse_dfs().
	
	Input:
		* args - parsed command line arguments
	'''
	
	sys.stdout.write('Parsing expression CSV file...\n')
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
	#spit it back out to function
	return expr_df

def filter_dfs(expr_df, deg_df):
	'''
	Take the expr_df and deg_df and trim them down to genes DE in at least one condition, and unify condition span.
	
	Input:
		* expr_df - PANDAS data frame of expression
		* deg_df - PANDAS data frame of DEG status
	'''
	
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
	
	return (expr_df, deg_df)

def sanity_checks(expr_df, deg_df, args):
	'''
	Make sure everything is okay with the input at this stage and nothing will trip up.
	'''
	
	#need at least 2 conditions (and preferably more than 2 genes too)
	if 0 in deg_df.shape or 1 in deg_df.shape:
		sys.stderr.write('ERROR: Insufficient number of conditions and/or genes\n')
		sys.exit(1)
	
	#need some genes to mine as well
	if all(np.sum(deg_df.values,axis=1)<2):
		sys.stderr.write('ERROR: No genes differentially expressed across two or more conditions\n')
		sys.exit(1)
	
	#need unique identifiers
	if (len(set(expr_df.index)) != len(list(expr_df.index))) or (len(set(deg_df.index)) != len(list(deg_df.index))):
		sys.stderr.write('ERROR: Non-unique identifiers present\n')
		sys.exit(1)
	
	#also would be nice that if we have the thresholding sizes, they actually match in number to condition span
	if args.thresh:
		if len(args.thresh) != (len(deg_df.columns)-1):
			sys.stderr.write('ERROR: Incompatible number of thresholding set sizes and actual number of data sets\n')
			sys.exit(1)

def parse_dfs(args):
	'''
	Take all things command line and data and return nice, handy, PANDAS parsed versions of gene expression, DEGs and other arguments.
	
	Input:
		* args - parsed command line arguments
	'''

	#parse the expression CSV
	expr_df = parse_expr_df(args)

	#parse the DEG CSV
	if args.deg:
		sys.stdout.write('Parsing DEG CSV file...\n')
		deg_df = pd.read_csv(args.deg, index_col=0, header=0)
		#once again, caps gene names on sight
		deg_df.index = map(lambda x:x.upper(), deg_df.index.tolist())
		deg_df.columns = map(lambda x:x.title(), deg_df.columns.tolist())
		#filter contents
		(expr_df, deg_df) = filter_dfs(expr_df, deg_df)
	else:
		#just make a mock one if there is none, makes it easier for the script later on
		#no need for filtering as obviously all is well
		sys.stdout.write('Creating mock DEG dataframe...\n')
		expgenes = expr_df.index.tolist()
		expconds = list(set(expr_df.columns.get_level_values('Treatment')))
		deg_df = pd.DataFrame(np.ones((len(expgenes),len(expconds))), index=expgenes, columns=expconds)

	#time for a few sanity checks.
	sanity_checks(expr_df, deg_df, args)

	#time for data standardisation
	if args.stand:
		sys.stdout.write('Standardising data...\n')
		for cond in list(set(expr_df.columns.get_level_values('Treatment'))):
			vals = expr_df.loc[:,cond].values
			expr_df.loc[:,cond] = (vals - np.mean(vals,axis=1)[:,None]) / np.std(vals,axis=1)[:,None]
	
	#I think we have everything this function can do. spit it out for further consumption
	return (expr_df, deg_df)

def wigwams_analysis_default(expr_df, deg_df, args):
	'''
	The bog standard Wigwams pipeline, as described in the paper: mining -> merging -> sweeping
	
	Input:
		* expr_df - PANDAS data frame of gene expression, as parsed by parse_dfs(). Multiindexing of columns - treatments and times.
		* deg_df - PANDAS data frame of DEG status, as parsed by parse_dfs().
		* args - command line arguments, as parsed by parse_dfs().
	'''

	#mining. the time consuming part
	#ww.mining(expr_df, deg_df, pool=args.pool, sets=args.sets, alpha = args.alpha, corrnet=args.corrnet, legacy=args.leg, job=args.job)
	#merging
	#ww.merging(expr_df, overlap=args.merging_overlap, meancorr=args.meancorr, corrfilt=args.corrfilt, legacy=args.leg, job=args.job)
	#sweeping
	#ww.sweeping(overlap=args.sweeping_overlap, job=args.job)
	#thresholding if applicable
	which_file = 'filtered'
	if args.thresh:
		ww.thresholding(sizes=args.thresh, job=args.job)
	else:
		which_file = 'swept'
	#exporting
	ww.export(expr_df, deg_df, which_file=which_file, annot_file=args.annot, hyper=args.hyper, stand=args.stand, job=args.job)

def run_single_condspan(expr_df, deg_df, args, swov):
	'''
	Helper function that runs a single condspan's worth of redundancy removal.
	'''
	
	mergefile = 'raw'
	sweepfile = 'merged'
	if args.thresh:
		sweepfile = 'filtered'
	for condspan in np.arange(len(deg_df.columns),1,-1):
		ww.merging(expr_df, overlap=args.merging_overlap, meancorr=args.meancorr, corrfilt=args.corrfilt, condspan=condspan, which_file=mergefile, legacy=args.leg, job=args.job)
		#need this mergefile trick so that we can call it on raw first and then on swept onwards
		mergefile = 'swept'
		if args.thresh:
			ww.thresholding(sizes=args.thresh, condspan=condspan, which_file='merged', job=args.job)
		#no point sweeping if we've only got two conditions' worth of stuff anymore
		if condspan > 2:
			ww.sweeping(overlap=swov, condspan=condspan, which_file=sweepfile, job=args.job)
	#assess how well we're doing
	#potential output is the current sweepfile
	modules = ww.read_modules(os.path.normcase(args.job+'/'+sweepfile+'_modules.tsv'))
	indrat = ww.module_size_check(modules)
	return indrat

def wigwams_analysis_new(expr_df, deg_df, args):
	'''
	The modified Wigwams pipeline, aiming to reduce redundancy across inter-condition modules and produce smaller modules across small condition spans.
	
	Pseudocode: mining -> for (i in number_of_conditions : -1 : 2): {merging(i) -> thresholding(i) -> sweeping(i)}
	
	Input:
		* expr_df - PANDAS data frame of gene expression, as parsed by parse_dfs(). Multiindexing of columns - treatments and times.
		* deg_df - PANDAS data frame of DEG status, as parsed by parse_dfs().
		* args - command line arguments, as parsed by parse_dfs().
	'''
	
	#mining. the time consuming part
	ww.mining(expr_df, deg_df, pool=args.pool, sets=args.sets, alpha = args.alpha, corrnet=args.corrnet, legacy=args.leg, job=args.job)
	#now, perform merging -> thresholding -> sweeping in succession for decreasing condition spans
	#loop over different "building block" sweeping stringencies
	swov_best = 0
	ratio_best = 0
	for swov in np.arange(args.sweeping_overlap,args.sweep_low-0.01,-0.05):
		indrat = run_single_condspan(expr_df, deg_df, args, swov)
		if indrat > ratio_best:
			ratio_best = indrat
			swov_best = swov
			sys.stdout.write('New least redundant output (with sweeping overlap '+str(swov)+')\n')
	#our export swoops in to pick up where the sweeping on two would have happened, so export the current sweepfile
	sys.stdout.write('Least redundant output obtained for sweeping overlap ratio of '+str(swov_best)+'.\n')
	#if the best values were obtained for the last one we checked, easy enough to just keep it
	if swov != swov_best:
		sys.stdout.write('Preparing output.\n')
		run_single_condspan(expr_df, deg_df, args, swov_best)
	sweepfile = 'merged'
	if args.thresh:
		sweepfile = 'filtered'
	ww.export(expr_df, deg_df, which_file=sweepfile, annot_file=args.annot, hyper=args.hyper, stand=args.stand, job=args.job)

def main():
	#parse the command line stuff
	args = parse_args()
	#parse the inputs
	(expr_df, deg_df) = parse_dfs(args)
	if not args.nro:
		#default analysis pipeline. split off into function to experiment with workflows if desired
		wigwams_analysis_default(expr_df, deg_df, args)
	else:
		#tweaked experimental pipeline
		wigwams_analysis_new(expr_df, deg_df, args)

#this is a pretty cool Sam trick
#if I put the whole code in main() then I can access it from something like iPython
#and at the same time have more functions in here if I so desire
#whilst this little bit here at the bottom makes it so that command line fires main()
if __name__ == "__main__":
	main()