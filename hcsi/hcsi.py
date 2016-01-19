import numpy as np
import scipy as sp
import multiprocessing as mp
import pandas as pd
import argparse
import sys
import warnings

import csi
import AbstractGibbs as ag
from main import parse_gp_hyperparam_priors

def loadData(fd):
	# read in data and make sure headers are correct
	inp = pd.read_csv(fd,dtype=str,index_col=0,header=[0,1])
	inp.columns = pd.MultiIndex.from_tuples([(a,float(b)) for a,b in inp.columns],names=inp.columns.names)
	# convert to floating point values
	return inp.astype(float)

def parse_args():
	#the obligatory argument parsing
	parser = argparse.ArgumentParser()
	parser.add_argument('--Input', dest='input', type=argparse.FileType('r'), required=True, help='The CSV file featuring your expression data. First column for gene names, first row for condition names repeated for each condition time point, second row for condition time point information.')
	parser.add_argument('--Depth', dest='depth', type=int, default=2, help='CSI parental set depth truncation to maintain computational tractability. Default: 2')
	parser.add_argument('--Prior', dest='gpprior', type=str, default='10,0.1', help='CSI Gaussian Process prior, provided as \'shape,scale\' for a gamma distribution or \'uniform\' for a uniform distribution. Default: \'10,0.1\'')
	parser.add_argument('--Genes', dest='genes', default=None, nargs='+', help='Child gene set to evaluate, if you wish to only run hCSI on a subset of the available gene space. Provide as space delimited names matching the CSV file. Default: None (analyse the whole dataset)')
	parser.add_argument('--Pool', dest='pool', type=int, default=1, help='Number of threads to open up for parallelising hCSI on a per-gene basis. Default: 1 (no parallelising)')
	parser.add_argument('--Samples', dest='samples', type=int, default=3000, help='Number of Gibbs updates to perform within hCSI. Default: 3000')
	parser.add_argument('--BurnIn', dest='burnin', type=int, default=500, help='Number of initial Gibbs updates to discard as burn-in. Default: 500')
	args = parser.parse_args()
	return args

def hamming(p1, p2):
	#filter out to what we actually want from the parents
	p1_set = set(p1[0])
	p2_set = set(p2[0])
	#the hamming distance is the total of elements that differ between the parent combos
	ulen = len(p1_set.union(p2_set))
	ilen = len(p1_set.intersection(p2_set))
	#so we return the total number that only show up in one of the two
	return ulen-ilen

def ismember(a, b):
    #returns True for the positions of a that are in b and false otherwise
    bind = {}
    for i, elt in enumerate(b):
        if elt not in bind:
            bind[elt] = True
    return np.asarray([bind.get(itm, False) for itm in a])

class GibbsHCSI(ag.AbstractGibbs):
    def __init__(self, rvObjList, rvHyperNetwork):
    #add a separate hypernetwork rv
        self.rvl = rvObjList
        self.indexList = range(len(self.rvl))
        self.sampledValues = []
        for i in self.indexList:
            self.sampledValues.append([])
        self.rvHyperNetwork = rvHyperNetwork
        self.sampledValuesHyperNetwork = []
    
    def gibbsUpdate(self):
    #expand the thing with hypernetwork sampling and parameter MCMC
        #run through all the variables
        valList = []
        betaList = []
        pset = self.rvHyperNetwork.getValRange()
        for i in self.indexList:
            #for each variable, build a distribution accross all possible values for that variable
            #based on current values of other variable 
            distribution = self.rvl[i].getConditionalDistribution(self.rvHyperNetwork.getCurrentValue())
            #sample a random value based on the distribution
            valInd = np.random.choice(np.arange(len(pset)),size=1,p=distribution)
            self.rvl[i].setCurrentValue(pset[valInd[0]])
            valList.append(pset[valInd[0]])
            betaList.append(self.rvl[i].beta)
        #sample the hypernetwork
        distribution = self.rvHyperNetwork.getConditionalDistribution(valList,betaList)
        valInd = np.random.choice(np.arange(len(pset)),size=1,p=distribution)
        self.rvHyperNetwork.setCurrentValue(pset[valInd[0]])
        '''TODO: theta and beta MCMC step'''

    def sample(self, repeats=1):
    #add hypernetwork sample storage
        for i in range(repeats):
            self.gibbsUpdate();
            for j in self.indexList:
                self.sampledValues[j].append(self.rvl[j].getCurrentValue())
            self.sampledValuesHyperNetwork.append(self.rvHyperNetwork.getCurrentValue())

class RandomVariableCondition(ag.RandomVariable):
    def __init__(self, csidata, cond, gene, gpprior, depth):
    #override init to include all sorts of CSI stuffs
        #extract the data frame columns that actually have our condition's data
        #(without killing the fine balance of the tuple indexing)
        colNames = np.asarray([x[0] for x in csidata.columns.values])
        numtemp = np.arange(len(colNames))
        inp = csidata.iloc[:,numtemp[colNames==cond]]
        #some processing stuff borrowed from the CSI code
        assert (inp.columns.is_monotonic_increasing)
        self.cc = csi.Csi(inp)
        self.em = self.cc.getEm()
        if gpprior:
            self.em.set_priors(gpprior[0], gpprior[1])
        #prepare the EM object
        self.em.setup(self.cc.allParents(gene,depth))
        '''TODO: sample beta'''
        self.beta = 1
        '''TODO: sample beta'''
        #going rogue here - starting with empty network
        self.currentValue = self.em.pset[0];
        self.valRange = self.em.pset;
        self.distribution = list(range(len(self.valRange)));
    
    def getConditionalDistribution(self, hyperparent):
    #override with actual computation
        i = 0
        logliks = self.em.logliks()
        #turn logliks into actual numbers, scaled to the largest one
        logliks = np.exp(logliks-np.max(logliks))
        #for each possible value of the RV
        for v in self.getValRange():
            #skipping P(theta) and all the Z's as they're the same for each parent combo
            p =  logliks[i] * np.exp((-1)*self.beta*hamming(v,hyperparent))
            #set the distribution
            self.distribution[i] = p
            i = i + 1
        self.normalizeDistribution()
        return self.distribution
    
    def normalizeDistribution(self):
    #we actually need the normalisation
        self.distribution = self.distribution/np.sum(self.distribution)

class RandomVariableHyperNetwork(ag.RandomVariable):
    def __init__(self, csidata, gene, depth):
    #override init to include all sorts of CSI stuffs
	    cc = csi.Csi(csidata)
	    self.valRange = cc.allParents(gene,depth);
	    self.currentValue = self.valRange[0];
	    self.distribution = list(range(len(self.valRange)));
    
    def getConditionalDistribution(self, condParents, betas):
    #override with actual computation
        i = 0
        for v in self.getValRange():
            p = 1
            for (j, cond) in enumerate(condParents):
                p = p * np.exp((-1)*betas[j]*hamming(v,cond))
            self.distribution[i] = p
            i = i + 1
        self.normalizeDistribution()
        return self.distribution

    def normalizeDistribution(self):
    #we actually need the normalisation
        self.distribution = self.distribution/np.sum(self.distribution)

def runGibbs(gene, inp, gpprior, args):
	hnrv = RandomVariableHyperNetwork(inp,gene,args.depth)
	crv = []
	conditions = np.unique([x[0] for x in inp.columns.values])
	for cond in conditions:
		crv.append(RandomVariableCondition(inp,cond,gene,gpprior,args.depth))
	gibbs = GibbsHCSI(crv,hnrv)
	gibbs.sample(args.samples)
	#ditch the burn-in
	for i in np.arange(len(gibbs.sampledValues)):
		gibbs.sampledValues[i] = gibbs.sampledValues[i][args.burnin:]
	gibbs.sampledValuesHyperNetwork = gibbs.sampledValuesHyperNetwork[args.burnin:]
	#parse rows of the magical marginal matrices
	out = np.zeros((len(conditions)+1,len(args.genes)))
	numtemp = np.arange(len(args.genes))
	for i in np.arange(len(conditions)):
		for j in np.arange(len(gibbs.sampledValues[i])):
			mask = ismember(args.genes, list(gibbs.sampledValues[i][j][0]))
			out[i,numtemp[mask]] += 1
		out[i,:] /= len(gibbs.sampledValues[i])
	#repeat for the hypernetwork
	for i in np.arange(len(gibbs.sampledValuesHyperNetwork)):
		mask = ismember(args.genes, list(gibbs.sampledValuesHyperNetwork[i][0]))
		out[-1,numtemp[mask]] += 1
	out[-1,:] /= len(gibbs.sampledValuesHyperNetwork)
	return (out, gene)

def _pool_init(inp_, gpprior_, args_):
	global inp, gpprior, args
	inp = inp_
	gpprior = gpprior_
	args = args_

def pool_runGibbs(gene):
	return runGibbs(gene,inp, gpprior, args)

def main():
	#read arguments
	args = parse_args()
	#parse data
	inp = loadData(args.input)
	conditions = np.unique([x[0] for x in inp.columns.values])
	#assess indegree sanity
	if args.depth < 1:
		sys.stderr.write("Error: truncation depth must be greater than or equal to one")
		sys.exit(1)
	#parse GP priors
	if args.gpprior is None or args.gpprior == 'uniform':
		gpprior = None
	else:
		try:
			gpprior = parse_gp_hyperparam_priors(args.gpprior)
		except ValueError(s):
			sys.stderr.write("Error: "+s)
			sys.exit(1)
	#parse gene list
	genes = args.genes
	if genes is None:
		genes = list(inp.index)
		#not pretty, but passes it on to the parallelised workers with ease
		args.genes = genes
	else:
		missing = np.setdiff1d(genes, inp.index)
		if len(missing) > 0:
			sys.stderr.write("Error: The following genes were not found: {missing}\n".format(missing=', '.join(missing)))
			sys.exit(1)
	#prepare output things
	output = []
	for i in range(len(conditions)+1):
		template = np.ndarray((len(genes)+1,len(genes)+1),dtype=object)
		template[0,0] = 'ID'
		template[1:,0] = np.asarray(genes)
		template[0,1:] = np.asarray(genes)
		output.append(template)
	numtemp = 1+np.arange(len(genes))
	#"Parallelised or not, here I come!" - Gibbs sampler, 2016
	if args.pool > 1:
		p = mp.Pool(args.pool, _pool_init, (inp, gpprior, args))
		for (out, gene) in p.imap_unordered(pool_runGibbs, genes):
			#wrap it into a one-element list so that ismember sees it whole
			mask = ismember(genes,[gene])
			ind = numtemp[mask][0]
			#sneaking in the hypernetwork as the last "condition"
			for i in np.arange(len(conditions)+1):
				output[i][ind,1:] = np.asarray([str(x) for x in out[i,:]])
	else:
		for gene in genes:
			#we don't need to re-catch the gene as we have it already
			(out, blaa) = runGibbs(gene, inp, gpprior, args)
			#wrap it into a one-element list so that ismember sees it whole
			mask = ismember(genes,[gene])
			ind = numtemp[mask][0]
			#sneaking in the hypernetwork as the last "condition"
			for i in np.arange(len(conditions)+1):
				output[i][ind,1:] = np.asarray([str(x) for x in out[i,:]])
	#spit the things out
	conditions = list(conditions)
	conditions.append('hypernetwork')
	for i in np.arange(len(conditions)):
		writer = open('hcsi_'+conditions[i]+'.csv','w')
		for j in np.arange(output[i].shape[0]):
			writer.write(','.join(output[i][j,:])+'\n')
		writer.close()

if __name__ == '__main__':
    main()