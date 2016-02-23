import random
from datetime import datetime


import numpy;
import pdb;

class AbstractGibbs(object):
    #rvObjList is a list of random variables object, each object has
    #a function to compute its conditional probability given values for other
    #variable in the list
    #sampledValues is a list of lists where each 'inner' list is the value of a varialbe sampled by the sampler at each update step
    def __init__(self, rvObjList):
        self.rvl = rvObjList;
        self.indexList = range(len(self.rvl));
        self.sampledValues = [];
        for i in self.indexList:
            self.sampledValues.append([]);
    
    #update step
    def gibbsUpdate(self, c, index, db = False):
        for i in self.indexList:
            #for each variable, build a distribution accross all possible values for that variable
            #based on current values of other variable
            distribution, scores = self.rvl[i].getConditionalDistribution(self.rvl[:i]+self.rvl[i+1:], db);
            #sample a random value based on the distribution
            #val = mysample(self.rvl[i].getValRange(),size=1,p=distribution);
            #self.rvl[i].setCurrentValue(val[0]);
            counter = 0;
            tot = 0;
            val = None;
            for d in distribution:
                if tot + d > c[index]:
                    temp = self.rvl[i].getValRange();
                    val = temp[counter];
                    break;
                else:
                    tot = tot + d;
                    counter = counter + 1;
            self.rvl[i].setCurrentValue(val);
            index += 1;
            if index == len(c):
                c = numpy.random.random(1000);
                index = 0;
        return index,c;
            

    #run Gibbs udpate with with "repeats" number of cycles        
    def sample(self, repeats=1):
        rc = numpy.random.random(1000);
        index = 0;
        for i in range(repeats):
            if (i + 1) % 100 == 0:
                print("Done "+str(i+1));
            index,rc = self.gibbsUpdate(rc,index,True);        
            for j in self.indexList:
                self.sampledValues[j].append(self.rvl[j].getCurrentValue());
                


#class RandomVariable has 4 important methods:
#1) conditionalProbability
#2/3) get/setCurrentValue
#4) getConditionalDistribution
class RandomVariable(object):
        #currentValue is the value that the variable is taking on
        #valRange is a range of all possible values that this variable can take
        #normalize should probably be set to true, to make sure that distribution is always between [0,1]
        #distribution is a list of weight corresponding to each value in valRange
    def __init__(self,val = None, valRange = None, normalize = False):
        self.currentValue = val;
        self.valRange = valRange;
        self.normalize = normalize;
        self.distribution = [];
    #getter for valRange
    def getValRange(self):
        return self.valRange;
    #getter for currentValue
    def getCurrentValue(self):
        return self.currentValue;
        
    #set currentValue and do any updates that might come afterwards
    def setCurrentValue(self, val):
        self.currentValue = val;
        self.updateAfterSettingValue(val);
        
    #to be implemented by anything extending these clases
    def updateAfterSettingValue(self, val):
        pass
    #compute the conditional probability of current RV given other RVs (and their values)        
    def conditionalProbability(self,value,rvList):
        pass
    
    #compute distribution
    def getConditionalDistribution(self,rvList, db):
        i = 0;
        self.distribution = [];
        #for each possible value of the RV
        for v in self.getValRange():
            #compute the conditional probability that this RV takes this value given other RVs (and their values)
            p =  self.conditionalProbability(v,rvList) ;
            #set the distribution
            self.distribution.append(p);
            i = i + 1;        
        debugscore = self.normalizeDistribution(db);                
        return self.distribution, debugscore;
    
    def normalizeDistribution(self, db):
        pass;
