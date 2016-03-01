##ROC
#helper class for the evaluation of ROC curves and area under ROC

import scipy as S
import numpy as N
import pdb


def roc(labels, predictions):
    """roc - calculate receiver operator curve
    labels: true labels (>0 : True, else False)
    predictions: the ranking generated from whatever predictor is used"""
    #1. convert to arrays
    labels = S.array(labels).reshape([-1])
    predictions = S.array(predictions).reshape([-1])

    #threshold
    t = labels>0
    
    #sort predictions in desceninding order
    #get order implied by predictor (descending)
    Ix = S.argsort(predictions)[::-1]
    #reorder truth
    t = t[Ix]

    #compute true positiive and false positive rates
    tp = S.double(N.cumsum(t))/t.sum()
    fp = S.double(N.cumsum(~t))/(~t).sum()

    #add end points
    tp = S.concatenate(([0],tp,[1]))
    fp = S.concatenate(([0],fp,[1]))

    return [tp,fp]

def auroc(labels=None,predictions=None,tp=None,fp=None):
    """auroc - calculate area under the curve from a given fp/rp plot"""

    if labels is not None:
        [tp,fp] = roc(labels,predictions)
    n = tp.size
    auc = 0.5*((fp[1:n]-fp[0:n-1]) * (tp[1:n]+tp[0:n-1])).sum()
    return auc
    pass


    

    
