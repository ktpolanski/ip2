"""Expression TwoSample test with GPs"""

#path for pygp, stats
import sys
sys.path.append('./../')
import pygp.covar.sederiv as sederiv
import pygp.covar.noiseCF as noiseCF
import pygp.covar.sumCF as sumCF
import pygp.gpr as GPR
import pygp.gprEP as GPREP
import pygp.EPLikelihood as EPL
import scipy as S
import pylab as PL
import pdb
import copy
import logging as LG
#log gamma priors for GP hyperparametersf
from mlib.stats.lnpriors import *


class TSData(object):
    """TSData
    - class for two sample test data
    """

    def __init__(self):
        pass


    def getExpr(self,name,series):
        """getExpr(name,series)
         - name: name of expresion set(gene)
         - series: series, i.e. normally either 0/1"""
        print "abstract class, implement me"
        return None

        


class GPTwoSampleMLII(object):
    __slots__=["covar","logtheta0","priors","maxiter","opt","Smean","likelihood","robust","gpr_0","gpr_1","gpr_s","gpr_join","likelihood"]
    """GP Two sample test where hyper parameters are optimized via ML type II"""

    def __init__(self,covar=None, logtheta0=None, maxiter=20,priors = None, dimension=1,opt=True,Smean=True,robust=False):       
        """__init__(data,covar=None,logtheta0=None,maxiter=20,,priors=None,dimension=1)
        - data: TSdata object
        - covar: covariance function (instance)
        - logtheta0: start point for optimzation of covariance function
        - maxiter:  iterations for optimization
        - priors:   list with priors for each hyper parameter
        - dimension: data dimension(1)
        - opt: perform optimization?
        - Smean: subtract mean ?
        - robust: robust likelihood (false)
        """

        def getDefaultCovar():
            """create default covarince function"""
            x_dim = 1
            #se 
            se_cf = sederiv.SquaredExponentialCFnn(index = arange(0,x_dim))
            #noise
            noise_cf = noiseCF.NoiseCovariance()
            covar = sumCF.SumCovariance(covars=[se_cf,noise_cf])

            #prior
            priors = []
            #1. amplitude
            priors.append([lngammapdf,[1,1]])
            #2. lengthscale
            scale = 30.0
            Lmin  = 3.0   #it's all rescaled -5..5

            loc   = Lmin/scale
            priors.append([lngammapdf,[scale,loc]])
            #3. noise
            Lmin = 0.5
            scale = 1.0

            loc   = Lmin/scale
            priors.append([lngammapdf,[scale,loc]])
            #logtheta0
            #set to mean for all priors
            logtheta0 = S.array([S.log(p[1][0]*p[1][1]) for p in priors])
            return [covar,logtheta0,priors]

        def getRobustCovar():
            """covariance structure for robust regression"""
            [covar,logtheta0,priors] = getDefaultCovar()
            #reduce covar to SE covar
            covar = covar.covars[0]

            sigma = S.exp(logtheta0[-1])
            c = 0.9
            logthetaL = S.log([c,1-c,sigma,1E4])
            #include likelihood covariance function
            logtheta0  = S.concatenate((logtheta0[0:-1],logthetaL))
            return [covar,logtheta0,None]

        
        #default Covar and prior:
        self.maxiter = maxiter
        #opt flag
        self.opt = opt
        self.Smean = Smean

        #defualt values for logtheta,prior and covar
        self.robust = robust
        if robust:
            self.likelihood = EPL.MOGLikelihood(Ncomponents=2)
            [_covar,_logtheta0,_priors] = getRobustCovar()
        else:
            [_covar,_logtheta0,_priors] = getDefaultCovar()

        if covar is None:
            covar = _covar
        if logtheta0 is None:
            logtheta0 = _logtheta0
        if priors is None:
            priors = _priors

        self.covar  = covar
        self.priors = priors
        self.logtheta0 = logtheta0

        pass


    def preprocess(self,M0,M1,logtrafo=False,rescale=False):
        """preprocessing.
        -this function can do log transormation and rescaling of raw data"""
        if logtrafo:
            M0[1] = [S.log(M0[1][i]) for i in range(len(M0[1]))]
            M1[1] = [S.log(M1[1][i]) for i in range(len(M1[1]))]

        if rescale:
            Y0   = S.concatenate(M0[1])
            IY0  = ~S.isnan(Y0)
            S0   = Y0[IY0].std()
            Y1   = S.concatenate(M1[1])
            IY1  = ~S.isnan(Y1)
            S1   = Y1[IY1].std()
            #rescale
            M0[1] = [M0[1][i]/S0 for i in range(len(M0[1]))]
            M1[1] = [M1[1][i]/S1 for i in range(len(M1[1]))]
        return [M0,M1]



    def test(self,M0,M1,verbose=False,opt=None,logtrafo=False,rescale=False):
        """test for differential expression
        M0: dataset in condition 0
        M1: dataset in condition 1
        verbose: create verbose plot
        opt: optimise hyper(true)
        logtrafo: takes logs first (False)
        rescale: rescale to unit varaince (False)
        """

        if opt is None:
            opt = self.opt

        [M0,M1] = self.preprocess(M0,M1,logtrafo=logtrafo,rescale=rescale)
        #0. model: both in one bucket:

        ##changed this here:
        #xjoin = S.concatenate((M0[0],M1[0]),axis=1)
        #yjoin = S.concatenate((M0[1],M1[1]),axis=1)
        #Mjoin = [xjoin,yjoin]

        #this version is compatible with lists:
        Mjoin = self.Mconcatenate(M0,M1)
        MJ    = Mjoin
        
        gpr_join = self.getGPR(Mjoin,opt=opt)

        #1. both separately:
        gpr_0    = self.getGPR(M0,opt=False)
        gpr_1    = self.getGPR(M1,opt=False)
        gprs     = GPR.GroupGP([gpr_0,gpr_1])

        if opt:
            
            Ifilter = S.ones_like(gpr_0.logtheta)
#            Ifilter[-1] = 0
            logtheta=GPR.optHyper(gprs,gpr_0.logtheta,Ifilter=Ifilter,priors=gpr_0.priors,maxiter=self.maxiter)
            gpr_0.logtheta=logtheta
            gpr_1.logtheta=logtheta

        LML_join = -1.0*gpr_join.lMl(gpr_join.logtheta,lml=True,dlml=False)
        LML_0    = -1.0*gpr_0.lMl(gpr_0.logtheta,lml=True,dlml=False)
        LML_1    = -1.0*gpr_1.lMl(gpr_1.logtheta,lml=True,dlml=False)
        
        #bayes factor
        ratio    = (LML_0+LML_1-LML_join)
        #store result structures in local object
        self.gpr_join = gpr_join
        self.gpr_s    = gprs
        self.gpr_0    = gpr_0
        self.gpr_1    = gpr_1

        if verbose:
            LG.debug('theta0: %s' % str(S.exp(gpr_0.logtheta)))
            LG.debug('theta1: %s' % str(S.exp(gpr_1.logtheta)))
            X = S.linspace(M0[0][0].min(),M0[0][0].max(),100)
            
            #plot
            PL.clf()
            PL.hold(True)

            LG.debug("theta(0):" + str(S.exp(gpr_0.logtheta)))
            LG.debug("theta(1):" + str(S.exp(gpr_1.logtheta)))
            LG.debug("theta: " + str(S.exp(gpr_join.logtheta)))

            #plot the GP erorrbars and means first:
            self.plotGPpredict(gpr_0,M0,X,{'alpha':0.1,'facecolor':'r'},{'linewidth':2,'color':'r'})
            self.plotGPpredict(gpr_1,M0,X,{'alpha':0.1,'facecolor':'g'},{'linewidth':2,'color':'g'})
            self.plotGPpredict(gpr_join,M0,X,{'alpha':0.1,'facecolor':'b'},{'linewidth':2,'color':'b'})
            for rep in range(len(M0[0])):
                PL.plot(M0[0][rep],M0[1][rep],'r.--')
            for rep in range(len(M1[0])):
                PL.plot(M1[0][rep],M1[1][rep],'g.--')
            PL.xlim((X.min(),X.max()))
            PL.title('%s: %.4f' % ('',ratio))
            PL.xlabel('Time/hr')
            PL.ylabel('Log expression level')
#            Ymax = MJ[1].max()
#            Ymin = MJ[1].min()
#            DY   = Ymax-Ymin
#            PL.ylim([Ymin-0.1*DY,Ymax+0.1*DY])
        return ratio


    def plotData(self,M0,M1,legend_loc='upper left'):
        """plotting function to plot underlying raw data"""
        P0 = []
        P1 = []
        for rep in range(len(M0[0])):
            P0.append(PL.plot(M0[0][rep],M0[1][rep],'r.--'))
        for rep in range(len(M1[0])):
            P1.append(PL.plot(M1[0][rep],M1[1][rep],'g.--'))
        X = S.linspace(M0[0][0].min(),M0[0][0].max(),100)
        PL.xlim((X.min(),X.max()))
        PL.xlabel('Time/h')
        PL.ylabel('Log expression level')
        #create legend
        PL.legend((P0[0],P1[1]),('Treatment','Control'),loc=legend_loc)
        

    def plotGPpredict(self,GP,M,X,format_fill,format_line):
        """plotGPpredict(GP,X,format_fill,format_line)
        - GP regression and plotting of std as +/2 2std and mean prediction"""
        #vectorized version of X for G if needed
        if len(X.shape)<2:
            Xv = X.reshape(X.size,1)
        else:
            Xv = X
        #regression:
        [p_mean,p_std] = self.regress(GP,M,X=Xv)
        #plot std
        Xp = concatenate((X,X[::-1]))
        Yp = concatenate(((p_mean+2*p_std),(p_mean-2*p_std)[::-1]))
        PL.fill(Xp,Yp,**format_fill)
        PL.plot(X,p_mean,**format_line)

    def Mconcatenate(self,M0,M1):
        """concatenate to Ms"""
        M_ = [[],[]]
        M_[0].extend(M0[0])
        M_[0].extend(M1[0])
        M_[1].extend(M0[1])
        M_[1].extend(M1[1])
        return M_
        

    def M2GPxy(self,M):
        """M2GPxy(M)
        - conversion function returning x/y in arrays suitable for the gpr module
        - the standard version of this helper function does not really do much, but see the cluster variant to find out what this is for"""
        x1 = M[0]
        y1 = M[1]
        #a version that can deal with lists
        x  = []
        y  = []
        for i in range(len(x1)):
            x.extend(x1[i])
            y.extend(y1[i])

        x  = S.array(x)
        y  = S.array(y)
#        x  = S.concatenate((x1),axis=0)
#        y  = S.concatenate((y1),axis=0)
        I_keep = (~S.isnan(x)) & (~S.isnan(y))
        x = x[I_keep]
        y = y[I_keep]
        x = x.reshape(x.size,1)
        return [x,y]

    def getGPR(self,M,opt=True):
        """create a GP object, optimize hyperparameters for M[0]: x , M[1]: multiple time series y"""

        [x,y] = self.M2GPxy(M)

        logtheta = self.logtheta0

        if self.robust:
            gpr = GPREP.GPEP(covar=self.covar,Nep=3,likelihood=self.likelihood,Smean=True,x=x,y=y)
            pass
        else:
            gpr = GPR.GP(self.covar,Smean=True,x=x,y=y)
        if opt:
            Ifilter = S.ones_like(logtheta)
            LG.debug("opt")
            LG.debug('priors: %s' % str(self.priors))
            logtheta=GPR.optHyper(gpr,logtheta,Ifilter=Ifilter,priors=self.priors,maxiter=self.maxiter)

        gpr.logtheta = logtheta
        gpr.priors   = self.priors
        return gpr
        

    def regress(self,gpr,M,X=None):
        """regress one of the series and return predicted mean/std and raw data"""

        x  = S.concatenate((M[0]),axis=0)
        y  = S.concatenate((M[1]),axis=0)

        if X is None:
            X = S.linspace(x.min(),x.max(),100)
            X = X.reshape(X.size,1)
        [p_mean,p_std] = gpr.predict(gpr.logtheta,X)
        return [p_mean,p_std]
    







class GPTwoSampleClusterMLII(GPTwoSampleMLII):
    """same as GPTwoSampleMLII but using the cluster covariance function with one noise level per cluster"""

    def __init__(self,data,covar=None, logtheta0=None, maxiter=20,priors = None, dimension=1,clusters=4):
        #call super
        super(GPTwoSampleClusterMLII,self).__init__(data,covar=covar,logtheta0=logtheta0,maxiter=maxiter,priors=priors,dimension=dimension)
        #then implement the few changes - basically just to accommodate the changed covaraince
        #(additional parameters)

        self.covar = sederiv.SquaredExponentialClusterCF(dimension,clusters)
        if priors is None:
            priors = []
            if 1:
                #informative gamma priors for length scales, amplitude and noise:#
                #0. amplitude:
                priors.append([lngammapdf,[1,1]])
                #length scale
                for i in range(dimension):
                    priors.append([lngammapdf,[1,1]])
                #noise, one for each cluster
                for ic in range(clusters):
                    priors.append([lngammapdf,[1,1]])
            else:
                #zero prior
                for i in range(self.covar.getNparams()):
                    priors.append([lnzeropdf,[10,0.1]])


        self.priors = priors
        self.logtheta0 = logtheta0


    def M2GPxy(self,M):
        """M2GPxy(M)
        - conversion function returning x/y in arrays suitable for the gpr module
        - the standard version of this helper function does not really do much, but see the cluster variant to find out what this is for"""
        x1 = M[0]
        y1 = M[1]


        x1 = x1.reshape([x1.shape[0],x1.shape[1],1])
        #add extra dimension for the indicator which cluster
        x1 = S.concatenate((x1,x1),axis=2)
        x1[0,:,1] = 0
        x1[1,:,1] = 1
        x1[2,:,1] = 2
        x1[3,:,1] = 3
        #now just concatenate along the first dimension
        y  = S.concatenate((y1),axis=0)
        x  = S.concatenate((x1),axis=0)
        return [x,y]

    def getGPR(self,M,opt=True):
        """create a GP object, optimize hyperparameters for M[0]: x , M[1]: multiple time series y"""

        [x,y] = self.M2GPxy(M)

        #initialize GP
        logtheta = self.logtheta0
        priors = self.priors

        #adjust priors ?
        if 0:
            Yexp = 0.5
            scale = 20
            loc   = Yexp/scale
            priors[0][1] = [scale,loc]
            Lmin = (x[:,0].max(axis=0)-x[:,0].min(axis=0))/4
            scale = 30
            loc   = Lmin/scale

            priors[1][1] = [scale,loc]
            sigma = 1.0
            scale = 30

            sigma = 0.5
            scale = 3
            loc   = sigma/scale
            

        gpr = GPR.GP(self.covar,Smean=True,x=x,y=y)
        pydb.set_trace()
        if opt:
            #opt filter
            Ifilter = S.ones_like(logtheta)
            logtheta=GPR.optHyper(gpr,logtheta,Ifilter=Ifilter,priors=priors,maxiter=self.maxiter)
        #save optimised hyperparameters
        gpr.logtheta = logtheta
        return gpr


    def regress(self,gpr,M,X=None):
        """regress one of the series and return predicted mean/std and raw data"""

        [x,y] = self.M2GPxy(M)

        if X is None:
            X = S.linspace(x[:,0].min(),x[:,0].max(),100)
            X = X.reshape(X.size,1)
        X = S.concatenate((X,X),axis=1)
        #add this because of cluster (noise is ignored anyhow)
        X[:,1] = 0
        
        [p_mean,p_std] = gpr.predict(gpr.logtheta,X)
        return [p_mean,p_std]
        
        
