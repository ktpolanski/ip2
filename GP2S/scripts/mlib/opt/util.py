#Optimization tools
import numpy as N

def check_grad(f,df,p):
    rv =  []
    eps = 0.001
    for i in range(len(p)):
        dp = N.zeros(len(p),'double')
        dp[i]=eps
        rv.append((f(p+dp)-f(p))/eps)
    rv = N.array(rv)
    rv-=N.array(df(p))
    return rv
