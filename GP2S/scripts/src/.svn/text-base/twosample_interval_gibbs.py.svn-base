"""Two sample test on intervals
- this version uses propper gibbs resampling of the indicator variables
"""

#path for pygp, stats
import sys
sys.path.append('./..')

from pygp.covar import *
import pygp.gpr as GPR
import scipy as S
import pylab as PL
#import pydb
import copy as CP
#log gamma priors for GP hyperparametersf
from mlib.stats.lnpriors import *
from twosample import *
from hinton import *




class GPTwoSampleInterval(GPTwoSampleMLII):
    __slots__=["prior_Z"]


    def __init__(self,prior_Z=0.5,**kwargin):
        """__init__(data,covar=None,logtheta0=None,maxiter=20,,priors=None,dimension=1)
        - data: TSdata object
        - covar: covariance function (instance)
        - logtheta0: start point for optimzation of covariance function
        - maxiter:  iterations for optimization
        - priors:   list with priors for each hyper parameter
        - dimension: data dimension(1)
        - opt: perform optimization?
        - Smean: subtract mean ?
        - robust: use robust likelihood (False)
        """
        GPTwoSampleMLII.__init__(self,**kwargin)
        #prior for missing proportions
        self.prior_Z  = S.array([1-prior_Z,prior_Z])
        pass


    def test_interval(self,M0,M1,verbose=False,opt=None,Ngibbs_iterations=10):
        """test with interval sampling test
        M0: first expresison condition
        M1: second condition
        verbose: produce plots (False)
        opt: optimise hyperparameters(True)
        Ngibbs_iterations (10)
        """

        def rescaleInputs(X,gpr):
            """copy implementation of input rescaling
            - this version rescales in a dum way assuming everything is float
            """
            X_ = (X-gpr.minX)*gpr.scaleX - 5
            return X_
        

        def sampleIndicator(I,take_out=True):
            """sample all indicators that are true in the index vector I
            take_out: take indices out of the dataset first (yes)
            """
            PZ = S.zeros([2,I.sum()])

            if take_out:
                IS_ = IS & (~I)
                IJ_ = IJ & (~I)
            else:
                IS_ = IS
                IJ_ = IJ
            #set datasets
            self.gpr_0.setData(S.concatenate((M0R[0][:,IS_]),axis=0).reshape([-1,1]),S.concatenate((M0R[1][:,IS_]),axis=1),process=False)
            self.gpr_1.setData(S.concatenate((M1R[0][:,IS_]),axis=0).reshape([-1,1]),S.concatenate((M1R[1][:,IS_]),axis=1),process=False)
            self.gpr_join.setData(S.concatenate((MJR[0][:,IJ_]),axis=0).reshape([-1,1]),S.concatenate((MJR[1][:,IJ_]),axis=1),process=False)

            
            Yp0 = self.gpr_0.predict(self.gpr_0.logtheta,XT[I],mean=False)
            Yp1 = self.gpr_1.predict(self.gpr_0.logtheta,XT[I],mean=False)
            Ypj = self.gpr_join.predict(self.gpr_0.logtheta,XT[I],mean=False)
            #compare to hypothesis
            #calculate log likelihoods for every time step under different models

            #D0  = -0.5*(M0R[1][:,I]-Yp0[0])**2/Yp0[1] -0.5*S.log(Yp0[1])
            #D1  = -0.5*(M1R[1][:,I]-Yp1[0])**2/Yp1[1] -0.5*S.log(Yp1[1])
            #DJ  = -0.5*(MJR[1][:,I]-Ypj[0])**2/Ypj[1] -0.5*S.log(Ypj[1])
            
            #robust likelihood
            c =0.9
            D0  = c    * S.exp(-0.5*(M0R[1][:,I]-Yp0[0])**2/Yp0[1])*1/sqrt(2*S.pi*Yp0[1])
            D0 += (1-c)* S.exp(-0.5*(M0R[1][:,I]-Yp0[0])**2/1E8)*1/sqrt(2*pi*S.sqrt(1E8))
            D0 = S.log(D0)
            D1  = c    * S.exp(-0.5*(M1R[1][:,I]-Yp1[0])**2/Yp1[1])*1/sqrt(2*S.pi*Yp1[1])
            D1 += (1-c)* S.exp(-0.5*(M1R[1][:,I]-Yp1[0])**2/1E8)*1/sqrt(2*pi*S.sqrt(1E8))
            D1 = S.log(D1)
            DJ  = c    * S.exp(-0.5*(MJR[1][:,I]-Ypj[0])**2/Ypj[1])*1/sqrt(2*S.pi*Ypj[1])
            DJ += (1-c)* S.exp(-0.5*(MJR[1][:,I]-Ypj[0])**2/1E8)*1/sqrt(2*pi*S.sqrt(1E8))
            DJ = S.log(DJ)
            
            DS  = D0.sum(axis=0) + D1.sum(axis=0)
            DJ  = DJ.sum(axis=0)
            ES  = S.exp(DS)
            EJ  = S.exp(DJ)
            PZ[0,:] = self.prior_Z[0]*EJ
            PZ[1,:] = self.prior_Z[1]*ES
            
            PZ      /= PZ.sum(axis=0)
            #sample indicators
            Z       = S.rand(I.sum())<=PZ[1,:]
            if(IS_.sum()==1):
                Z = True
            if(IJ_.sum()==1):
                Z = False
            return [Z,PZ]
            pass
        pass

        #1. use the standard method to initialise the GP objects
        ratio = self.test(M0,M1,verbose=verbose,opt=opt)
        GP0  = CP.deepcopy(self.gpr_0)
        GP1  = CP.deepcopy(self.gpr_1)
        GPJ  = CP.deepcopy(self.gpr_join)
        
        #2. initialise gibbs samplign for time-local approximation
        #get set of unique time coordinates
        XT = M0[0][0].reshape([-1,1])
        #rescale all datasets inputs
        MJ = [S.concatenate((M0[0],M1[0]),axis=0),S.concatenate((M0[1],M1[1]),axis=0)]
        MJR = CP.deepcopy(MJ)
        M0R = CP.deepcopy(M0)
        M1R = CP.deepcopy(M1)
        #rescale and 0 mean
        MJR[0] = rescaleInputs(MJ[0],self.gpr_join)
        M0R[0] = rescaleInputs(M0[0],self.gpr_0)
        M1R[0] = rescaleInputs(M1[0],self.gpr_1)
        MJR[1]-= self.gpr_join.mean
        M0R[1]-= self.gpr_0.mean
        M1R[1]-= self.gpr_1.mean
        
        #index vector assigning each time point to either expert
        IS = S.ones([XT.shape[0]],dtype='bool')
        IJ = S.ones([XT.shape[0]],dtype='bool')
        
        #1. sample all indicators conditioned on current GP approximation
        Z_ = sampleIndicator(S.ones(XT.shape[0],dtype='bool'),take_out=False)
        Z_ = S.random.rand(IS.shape[0])>0.5
        PZ= S.zeros([2,IS.shape[0]])
        Z = S.zeros([Ngibbs_iterations,IS.shape[0]],dtype='bool')
        #update the datasets in the GP, attention: rescaling might cause trouble here..
        IS[~Z_] = False
        IJ[Z_]  = False
        #sample indicators one by one
        for n in range(Ngibbs_iterations):
#            for i in random.permutation(IS.shape[0]):
            for i in S.arange(IS.shape[0]):
                #resample a single indicator
                I = S.zeros([IS.shape[0]],dtype='bool')
                I[i] = True
                [Z_,PZ_]  = sampleIndicator(I)
                #save prob. and indicator value
                PZ[:,i] = S.squeeze(PZ_)
                Z[n,i]  = S.squeeze(Z_)
                #update indicators
                IS[i] = Z_
                IJ[i] = ~Z_
            LG.debug("Gibbs iteration:%d" % (n))
            LG.debug(PZ)
            LG.debug(IS)
        n_ = round(Ngibbs_iterations)/2
        Z  = Z[n_::]
        PZ[1,:] = Z.mean(axis=0)
        PZ[0,:] = 1-PZ[1,:]
        GP0 = self.gpr_0
        GP1 = self.gpr_1
        GPJ = self.gpr_join
        if verbose:
            PL.clf()
            #1. plot the gp predictions
            ax1=PL.axes([0.1,0.1,0.8,0.7])
            Xt_ = S.linspace(0,XT[:,0].max()+2,100)
            Xt  = Xt_.reshape([-1,1])
            self.plotGPpredict(GP0,M0,Xt,{'alpha':0.1,'facecolor':'r'},{'linewidth':2,'color':'r'})
            self.plotGPpredict(GP1,M0,Xt,{'alpha':0.1,'facecolor':'g'},{'linewidth':2,'color':'g'})
            self.plotGPpredict(GPJ,M0,Xt,{'alpha':0.1,'facecolor':'k'},{'linewidth':2,'color':'k'})
            PL.plot(M0[0].T,M0[1].T,'r.--')
            PL.plot(M1[0].T,M1[1].T,'g.--')

            
            PL.xlim([Xt.min(),Xt.max()])
            #remove last ytick to avoid overlap
            yticks = ax1.get_yticks()[0:-2]
            ax1.set_yticks(yticks)
            xlabel('Time/h')
            ylabel('Log expression level')
            Ymax = MJ[1].max()
            Ymin = MJ[1].min()
            DY   = Ymax-Ymin
            PL.ylim([Ymin-0.1*DY,Ymax+0.1*DY])

            #now plot hinton diagram with responsibilities on top
            ax2=PL.axes([0.1,0.715,0.8,0.2],sharex=ax1)
#           ax2=PL.axes([0.1,0.7,0.8,0.2])
            #PL.plot(XT[:,0],Z[1,:])
            #swap the order of Z for optical purposes
            Z_= S.ones_like(PZ)
            Z_[1,:] = PZ[0,:]
            Z_[0,:] = PZ[1,:]
            hinton(Z_,X=M0[0][0])
            ylabel('diff.')
            #hide axis labels
            setp( ax2.get_xticklabels(), visible=False)
            #font size
            #setp( ax1.get_xticklabels(), fontsize=tickSize)
            #setp( ax1.get_yticklabels(), fontsize=tickSize)
            #setp( ax2.get_xticklabels(), fontsize=tickSize)
            #PL.ylim(Ymin-0.1*DY,Ymax+0.1*DY)
        return [ratio,PZ]
        pass

    
        


    def test_interval_old(self,gene_name,verbose=True,opt=None):
        """test for differential expression with clustering model
        - returns a data structure which reflects the time of sepataoin (posterior over Z)
        """

        def updateGP():
            """update the GP datasets and re-evaluate the Ep approximate likelihood"""
            #0. update the noise level in accordance with the responsibilities
            for t in range(T):
                XS[:,t:R*T:T,-1] = 1/(Z[1,t]+1E-6)
                XJ[:,t:R*T:T,-1] = 1/(Z[0,t]+1E-6)

            GPS.setData(XS,Y)
            #here we joint the two conditions
            GPJ.setData(S.concatenate(XJ,axis=0),S.concatenate(Y,axis=0))
            #1. set the data to both processes

        M0 = self.data.getExpr(gene_name,0)
        M1 = self.data.getExpr(gene_name,1)
        MJ = [S.concatenate((M0[0],M1[0]),axis=0),S.concatenate((M0[1],M1[1]),axis=0)]


        C  = 2              #conditions
        R  = M0[0].shape[0] #repl.
        T  = M0[0].shape[1] #time
        D  = 2              #dim.


        #Responsibilities: components(2) x time                        
        Z  = 0.5*S.ones((2,T))

        #Data(X/Y): conditions x replicates x time x 2D
        X  = S.zeros((C,R*T,D))
        Y  = S.zeros((C,R*T))
        #unique times
        XT = S.ones((T,2))
        XT[:,0] = M0[0][0,:]

        [x0,y0] = self.M2GPxy(M0)
        [x1,y1] = self.M2GPxy(M1)
        
        X[0,:,0:2] = x0
        X[1,:,0:2] = x1
        Y[0,:]     = y0
        Y[1,:]     = y1
        #create indicator vector to identify unique time points

        
        #create one copy of the input per process as this is used for input dependen noise
        XS       = X.copy()
        XJ       = X.copy()


        #initialize the two GPs
        if self.logtheta0 is None:
            logtheta = self.covar.getDefaultParams()

        #the two indv. GPs
        GP0 = GPR.GP(self.covar,Smean=self.Smean,logtheta=self.logtheta0)
        GP1 = GPR.GP(self.covar,Smean=self.Smean,logtheta=self.logtheta0)
        #the group GP summarising the two indiv. processes
        GPS = GPR.GroupGP([GP0,GP1])
        #the joint process
        GPJ = GPR.GP(self.covar,Smean=self.Smean,logtheta=self.logtheta0)
        #update the GP
        updateGP()


        debug_plot = True

        for i in range(1):
            ###iterations###
            #1. get predictive distribution for both GPs
            ##debug
            #introduce the additional dimension to accom. the per obs. noise model
            #get prediction for all time points
            Yp0 = GP0.predict(GP0.logtheta,XT)
            Yp1 = GP1.predict(GP1.logtheta,XT)
            Ypj = GPJ.predict(GPJ.logtheta,XT)
            #considere residuals
            D0  = ((M0[1]-Yp0[0])**2 * (1/Yp0[1])).sum(axis=0)
            D1  = ((M1[1]-Yp1[0])**2 * (1/Yp1[1])).sum(axis=0)
            DJ  = ((MJ[1]-Ypj[0])**2 * (1/Ypj[1])).sum(axis=0)
            #the indiv. GP is the sum
            DS  = D0+D1
            #now use this to restimate Q(Z)
            ES  = S.exp(-DS)
            EJ  = S.exp(-DJ)
            #

            Z[0,:] =self.prior_Z[0]*EJ
            Z[1,:] =self.prior_Z[1]*ES
            Z     /=Z.sum(axis=0)
#             pydb.set_trace()
            updateGP()


        if verbose:
            PL.clf()
            labelSize = 15
            tickSize  = 12
            
            #1. plot the gp predictions
            ax1=PL.axes([0.1,0.1,0.8,0.7])
            Xt_ = S.linspace(0,XT[:,0].max()+2,100)
            Xt  = S.ones((Xt_.shape[0],2))
            Xt[:,0] = Xt_

            self.plotGPpredict(GP0,M0,Xt,{'alpha':0.1,'facecolor':'r'},{'linewidth':2,'color':'r'})
            self.plotGPpredict(GP1,M0,Xt,{'alpha':0.1,'facecolor':'g'},{'linewidth':2,'color':'g'})
            self.plotGPpredict(GPJ,M0,Xt,{'alpha':0.1,'facecolor':'b'},{'linewidth':2,'color':'b'})
            PL.plot(M0[0].T,M0[1].T,'r.--')
            PL.plot(M1[0].T,M1[1].T,'g.--')
            
            PL.xlim([Xt.min(),Xt.max()])
            #remove last ytick to avoid overlap
            yticks = ax1.get_yticks()[0:-2]
            ax1.set_yticks(yticks)
            xlabel('Time/h',size=labelSize)
            ylabel('Expression level',size=labelSize)

            #now plot hinton diagram with responsibilities on top
            ax2=PL.axes([0.1,0.715,0.8,0.2],sharex=ax1)
#            ax2=PL.axes([0.1,0.7,0.8,0.2])
            #PL.plot(XT[:,0],Z[1,:])
            #swap the order of Z for optical purposes
            Z_= S.ones_like(PZ)
            Z_[1,:] = Z[0,:]
            Z_[0,:] = Z[1,:]
            hinton(Z_,X=M0[0][0])
            ylabel('diff.')
            
            #hide axis labels
            setp( ax2.get_xticklabels(), visible=False)
            #font size
            setp( ax1.get_xticklabels(), fontsize=tickSize)
            setp( ax1.get_yticklabels(), fontsize=tickSize)
            setp( ax2.get_xticklabels(), fontsize=tickSize)
            #axes label
        return Z
        pass




