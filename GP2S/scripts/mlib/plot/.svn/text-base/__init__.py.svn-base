##IMAGE
#python module for matlab-compatible plotting:

import pylab
import pdb
import scipy.stats
import scipy as S

def imagesc(X,colormap=None):
    if(colormap==None):
        colormap = pylab.cm.gray
    pylab.imshow(X,aspect='auto',interpolation='nearest',cmap=colormap)
    
def savefig(fname):
    pylab.savefig("%s.eps" % (fname))
    pylab.savefig("%s.png" % (fname))    


def density_hist(Y,bins,X,*arg_in,**kwarg_in):
    """plot histogram and kernel density estimation in one plot
    Y: data to histogramify
    X: area to evaluate the kernel
    bins: bins argument for hist
    *argin, **kw_argin: parameters for density plot
    """
    Y = S.array(Y,dtype='float')
    Y_ = S.concatenate((Y,Y,Y))
#    Y_ = Y
    K = scipy.stats.gaussian_kde(Y_)
    #1. hist
    p0=pylab.hist(Y,bins)[2]
    #2. plot density
    Yk = K.evaluate(X)
    Yk*= Y.size
    p1=pylab.plot(X,Yk,*arg_in,**kwarg_in)
    return [p0,p1]
