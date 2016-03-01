# we don't have to do anything here.. but python needs this stupid file!

import numpy.random as R
import numpy as N


# implement mixture of gaussians

class MOG:
    "mexture of gaussians"
    def rvs(prob,params):
        '''rvs(prob,params): prob vector with prob of the components
                                            params: parameter of gaussians(cp numpy.random.multivaraite_normal)'''
        #first decide about which normal to use
        r = R.rand()
        i=0
        s=prob[0]
        while(r>s):
            i+=1
            s+=prob[i]
        #draw a sample from the multivariate normal
        #NOTE: multivariate-normal seem not work properly!
        return R.normal(params[i][0],params[i][1])
    rvs = staticmethod(rvs)
    def pdf(prob,params):
        "todo"
        pass

    
