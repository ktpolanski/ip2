#CLASS RV
#RandomVariable Class
#takes a distribution and it's parameters and generates transparently samples from it.. 
class RV:
    "Class representing a Random variable"
    T = None
    dist = None
    
    #:init(dist)
    def __init__(self,dist,T=[]):
        self.T = T
        self.dist = dist
        return
    #::init
    #:rvs
    def rvs(self,T=None,**kwds):
        "Generate a sample from the distribution with specified parameters"
        if(T==None):
            T  = self.T
        return float(self.dist.rvs(T[0],T[1]))
        #TODO: as soon as we pass on the kwds-filed we get a lot of crap as answer :(
    #:pdf
    def pdf(self,T=None):
        "Just return value odf PDF of the distribution with specified parameters"
        if(T==None):
                T  = self.T
        return self.dist.rvs(T[0],T[1])
    # ::sample
#::class RV


# implementation of MOG
