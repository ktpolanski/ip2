"""WarWick two sample test on bio data"""

import sys
sys.path.append('./../../')
sys.path.append('/home/os252/work/lib/python/mlib')

import pydb
import scipy as S
import mlib.io.csv as CSV
import cPickle
import os
from twosample import *
from twosample_interval import *
import pylab as PL
import pygp.gpr as GPR


#dump file for speedup:
dump_file   = '/home/os252/work/GP/GPTwoSample/data/warwick.dump'


class WarWickTSData(TSData):
    """WarWickTSData
    - data class for WarWick two sample experiment"""
    expr_matrix = None
    expr_header = None
    gene_names  = None

    def __init__(self,gene_file,expr_file):
        """__init__(gene_file,expr_file)
          - expr_file/gene_file from Warwick data"""

        #0. load genes and raw exp
        expr  = CSV.CSVFile(expr_file,delimiter=',',dtype='str',rowHeader=True,colHeader=False)
        genes = CSV.CSVFile(gene_file,delimiter=',',dtype='str',rowHeader=False,colHeader=False)
        self.expr_matrix = S.array(expr.getM(),dtype='double')
        self.expr_header = expr.getRowHeader()
        self.gene_names  = genes.getM()
        #set index for series 0/1 - there are 4 replicates each so we have quite a longish index list
        self.index = []
        self.index.append(S.zeros((4,24),dtype='int'))
        self.index.append(S.zeros((4,24),dtype='int'))
        for i in range(4):
            self.index[0][i,:] = 24*2*i     + S.arange(0,24)
            self.index[1][i,:] = 24*(2*i+1) + S.arange(0,24)

        #samling times:
        self.x = S.arange(2,50,2)
        pass


    def getExpr(self,name,series):
        """getExpr(name,series)
        -name: gene name
        -seris: series, i.e. normally either 0/1"""
        
        #0. get index of gene with name
        gene_index = S.where(self.gene_names==name)[0]
        #1. get data element
        expr   = self.expr_matrix[gene_index,:].T
        rv     = expr[self.index[series]].squeeze()
        #rv2    = self.expr_header[self.index[series]]
        #rv[3,:] = rv[3,:] +  5*S.random.randn(24)
        # create x array (multiple copies)
        x      = S.zeros_like(rv)
        for i in range(x.shape[0]):
            x[i,:] = self.x
        return [x,rv]
        
        



if __name__== '__main__':
    if(os.path.exists(dump_file)):
        data = cPickle.load(open(dump_file,"rb"))
    else:
        data = WarWickTSData(gene_file,expr_file)
        cPickle.dump(data,open(dump_file,'wb'),-1)

    
    GPR.DEBUG = 2

    gene_names = ['CATMA3A53880']

    logtheta0 = S.log([0.5,6,0.4])
    gptest = GPTwoSampleMLII(covar=None,logtheta0=logtheta0,maxiter=20)

    
