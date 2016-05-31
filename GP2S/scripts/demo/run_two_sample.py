#configure plotter to save plots
import matplotlib
matplotlib.use('Agg')
#append directories to path
import sys
sys.path.append('./..')
#debugger
import pdb
#IO libraries for reading csv files
import mlib.io.csv as CSV
#scientific python
import scipy as SP
import numpy as np
import pandas as pd
import argparse
import multiprocessing as mp
#pylab - matlab style plotting
import matplotlib.pylab as PL
import os
#log level control
import logging as LG
#two_sample: smooth model implements standard test and time dependent model. 
from src.twosample_interval_smooth import *

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='expr', type=argparse.FileType('r'), help='CSV file featuring the expression profiles of the genes across the control and the treatment. Gene names in first column. First two rows for headers - row one for condition name (e.g. control/treated), row two for time point.')
    parser.add_argument('--Pool', dest='pool', type=int, default=1, help='Number of processes to use when parallelising on a per-gene basis. Default: 1 (no parallelisation)')
    parser.add_argument('--TL', dest='tl', action='store_true', help='Flag. If provided, the time-local version of GP2S will be ran in an attempt to identify the time of first differential expression.')
    parser.add_argument('--Theta0', dest='theta0', type=str, default='0.5;1;0.4', help='The starting values for the hyperparameters for the Gaussian processes, delimited by commas. Default: 0.5,1,0.4')
    parser.add_argument('--Maxiter', dest='maxiter', type=int, default=20, help='Number of optimisation iterations. Default: 20')
    parser.add_argument('--ThetaZ', dest='thetaz', type=str, default='0.7;2;1E-5', help='TL-specific input. The values of the hyperparameters for the Gaussian process used to smooth out the Z, not optimised within the code. Default: 0.7,2,1E-5')
    parser.add_argument('--PriorZ', dest='priorz', type=float, default=0.3, help='TL-specific input. Prior belief in differential expression. Default: 0.3')
    parser.add_argument('--FixZ', dest='fixz', type=int, default=2, help='TL-specific input. The first N time points will be initially fixed to non-DE in the sampler. Default: 2')
    parser.add_argument('--ThreshZ', dest='threshz', type=float, default=0.5, help='TL-specific input. Z\'s greater than or equal to this value get interpreted as DE, and the first occurrence of Z with this threshold cleared will be reported as the ToFDE. Default: 0.5')
    parser.add_argument('--GibbsN', dest='ngibbs', type=int, default=30, help='TL-specific input, controlling the number of Gibbs samples. Default: 30')
    parser.add_argument('--Predpoints', dest='predpoints', type=int, default=20, help='TL-specific input. At how many (linearly spaced from the first to last time point) points should the GP be evaluated? Default: 20')
    args = parser.parse_args()
    
    args.theta0 = args.theta0.split(';')
    args.thetaz = args.thetaz.split(';')
    for i in range(len(args.theta0)):
        args.theta0[i] = np.log(float(args.theta0[i]))
    for i in range(len(args.thetaz)):
        args.thetaz[i] = np.log(float(args.thetaz[i]))
    args.theta0 = np.array(args.theta0)
    args.thetaz = np.array(args.thetaz)
    args.fixz = SP.array(np.arange(args.fixz))
    
    return(args)

def read_csv(args):
    #heavily ripped from Wigwams
    R = pd.read_csv(args.expr, index_col=0, header=0)
    #for some reason, this particular CSV setup scrambles the data with a double header
    #so import with a single header and make a tuple index for later by hand
    time = R.iloc[0][:]
    R = R[1:][:]
    R.index = map(lambda x:x.upper(), R.index.tolist())
    old_colnames = R.columns.values.tolist()
    new_colnames = []
    for i, lab in enumerate(old_colnames):
        lab = lab.split('.')
        if len(lab)==1:
            lab.append('0')
        new_colnames.append((lab[0].title(), float(time[i]), float(lab[1])))
    R.columns = pd.MultiIndex.from_tuples(new_colnames,names=['Treatment','Time','RepThing'])
    return(R)

def sanity_checks(R):
    #extract out condition names and assess that we have nice uniform time point distributions
    condnames = np.array([x[0] for x in list(R.columns)])
    tps = np.array([x[1] for x in list(R.columns)])
    conds = np.unique(condnames)
    #commence sanity checks
    if len(conds)!=2:
        sys.stderr.write('ERROR: More than two treatment specifications detected. Exiting\n')
        sys.exit(1)
    if np.sum(condnames==conds[0])!=np.sum(condnames==conds[1]):
        sys.stderr.write('ERROR: Unbalanced number of data points between the two treatments. Exiting\n')
        sys.exit(1)
    #okay, so if we made it this far we have the same number of data points and only two condition names
    tp_conds = []
    for cond in conds:
        cond_tps = tps[condnames==cond]
        #we should have the same exact number of reps per time point
        #so if we count up how many reps we have per time point, there should only be one unique value
        tpholder, inverse = np.unique(cond_tps, return_inverse=True)
        if len(np.unique(np.bincount(inverse)))!=1:
            sys.stderr.write('ERROR: Non-uniformity of time points for replicates detected in condition '+cond+'. Exiting\n')
            sys.exit(1)
        #well, if not, then we're fine and can store information
        tp_conds.append(tpholder)
        Nrepl = np.bincount(inverse)[0]
    #one last sanity check - are the time points the same?
    if not np.array_equal(tp_conds[0],tp_conds[1]):
        sys.stderr.write('ERROR: Different time points specified across the two treatments. Exiting\n')
        sys.exit(1)

def run_gp2s(gene,R,args):
    #create gptest object
    gptest = GPTwoSampleInterval(covar=None,logtheta0=args.theta0,maxiter=args.maxiter,prior_Z=args.priorz)
    #parse up the time structure
    Tc = np.unique([x[1] for x in list(R.columns)])
    Nrepl = R.shape[1]/(2*len(Tc))
    T  = SP.zeros([Nrepl,Tc.shape[0]])
    T[:,:] = Tc
    #get the condition names as they appear in R, with Y0 being the first and Y1 being the second
    conds = [R.columns[0][0], R.columns[-1][0]]
    #get data and reshape it
    Y0_hold = R.loc[gene][conds[0]].values
    Y1_hold = R.loc[gene][conds[1]].values
    Y0 = np.reshape(Y0_hold,(Nrepl,len(Tc)),order='C')
    Y1 = np.reshape(Y1_hold,(Nrepl,len(Tc)),order='C')
    #prepare GP2S structures
    M0 = [T,Y0]
    M1 = [T,Y1]
    if args.tl:
        #creates score and time local predictions
        Tpredict = SP.linspace(T.min(),T.max(),args.predpoints)[:,None]
        [score,Z] = gptest.test_interval(M0,M1,verbose=True,opt=True,Ngibbs_iterations=args.ngibbs,XPz=Tpredict,logthetaZ=args.thetaz,fix_Z=args.fixz) 
        Z1 = [str(x) for x in Z]
        str_out = gene+'\t'+str(score)
        tofde = np.argmax(Z>=args.threshz)
        #for reasons I don't get, it spits out 0 in two circumstances
        #when the condition is already satisfied at the start
        #and when it's not satisfied at all
        if tofde==0 and Z[0]<args.threshz:
            str_out = str_out+'\tNaN\n'
        else:
            #need to Tpredict[][] as it's a 2D'ish array
            str_out = str_out+'\t'+str(Tpredict[tofde][0])+'\n'
        str_out2=gene+'\t'+'\t'.join(Z1)+'\n'
    else:
        #only score
        score = gptest.test(M0,M1,verbose=True,opt=True)
        str_out = gene+'\t'+str(score)+'\n'
        str_out2=None
    #update figure window
    PL.draw()
    PL.savefig(os.path.join('./plots','%s.png'%gene))
    return(str_out,str_out2)

def _pool_init(R_, args_):
    global R, args
    R = R_
    args = args_

def run_gp2s_pool(gene):
    return run_gp2s(gene,R,args)

def main():
    #sponge up the command line arguments
    args = parse_args()

    #read CSV file
    R = read_csv(args)
    
    #sanity checks
    sanity_checks(R)
    
    #run GP2S proper
    writer = open('scores.txt','a')
    if args.tl:
        writer2 = open('Z.txt','a')
    if args.pool > 1:
        p = mp.Pool(args.pool, _pool_init, (R, args))
        for line, line2 in p.imap(run_gp2s_pool,list(R.index)):
            writer.write(line)
            if line2 is not None:
                writer2.write(line2)
    else:
        for gene in list(R.index):
            line, line2 = run_gp2s(gene,R,args)
            writer.write(line)
            if line2 is not None:
                writer2.write(line2)
    writer.close()
    if args.tl:
        writer2.close()

if __name__== '__main__':
    main()