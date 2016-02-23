from Network import Vertex as Vertex;
from AbstractGibbs import RandomVariable as RandomVariable;
from WLKernel import WLKernel as WLKernel;
import numpy as np;
import pdb;
import random;
class Gene (Vertex):
    #Each gene is a vertex in the GRN, labels mark the orthorlogues, MUST BE INTEGER!!!!!
    #gid is the id of each gene, could be string
    #speciesid is the id of each species, which is a number in 0,1,2,...,#species - 1
    #each gene is initialized as standing alone, not connected to any other genes
    def __init__(self, label, gid, speciesid):
        super(Gene,self).__init__(label,[],[]);
        self.geneId = gid;
        self.speciesId = speciesid;
    def getSpeciesId(self):
        return self.speciesId;
    def getGeneId(self):
        return self.geneId;
    def __str__(self):
        return 's'+str(self.speciesId)+'g'+str(self.geneId);
#look up results from CSI
class ParentWeightMap:
    def __init__(self,weightDictionary):
        self.weights = weightDictionary;
        
    def getWeight(self, currGene, geneList):
        key = str(currGene.getSpeciesId())+'#'+str(currGene.getGeneId());
        if len(geneList) > 0:
            paids = [pa.getGeneId() for pa in geneList];
            paids.sort();
            for pid in paids:
                key += '#'+str(pid);
        return self.weights[key];
    def isIn(self, currGene, geneList):
        key = str(currGene.getSpeciesId())+'#'+str(currGene.getGeneId());
        for pa in geneList:
            key += '#'+str(pa.getGeneId());
        return (key in self.weights);

class ParentSetRV(RandomVariable):
    def __init__(self, currGene, valRange, currNetwork, allNetworks, weightMap, temperature,h):
        super(ParentSetRV, self).__init__(currGene.getParents(), valRange, True);
        self.currNet = currNetwork;
        self.allNets = allNetworks;
        self.currGene = currGene;
        self.weightMap = weightMap;
        self.beta = temperature;
        self.h = h;
    def __str__(self):
        return self.currGene.__str__();    
        #update parental set of the currGene node in the graph (by adding/removing edges)
    def updateAfterSettingValue(self, val):
        self.currNet.updateParent(self.currGene, val)
    
    def conditionalProbability(self, value, rvList):
        pass
    
    def normalizeDistribution(self,db):
        scores = np.array([diff[0] for diff in self.distribution]);
        lookups = np.array([diff[1] for diff in self.distribution]);
        mmx = np.amax(scores);
        myE = mmx - scores;
        logprob = np.log(lookups);
        logDist = -1.0*self.beta*myE + logprob;
        mmx2 = np.amax(logDist);
        temp = np.exp(logDist - mmx2);
        sums = np.sum(temp);
        self.distribution = [t*1.0/sums for t in temp];
        return scores;

#Random variable to be used in 2nd framework, comparing graphs pair wise
class ParentSetRVFw2(ParentSetRV):
    def __init__(self,currGene,valRange,currNetwork,allNetworks,weightMap,temperature, h):
        super(ParentSetRVFw2,self).__init__(currGene,valRange, currNetwork, allNetworks, weightMap, temperature, h);
    def conditionalProbability(self, value, rvList):
        self.setCurrentValue(value);
        wlKernel = WLKernel(self.allNets);
        diffScore = wlKernel.computeWLKernelFramework2(self.h); 
        return [diffScore, self.weightMap.getWeight(self.currGene, self.getCurrentValue())];

#Random variable to be used in 1st framework, comparing each graph with hypernet    
class ParentSetRVFw1(ParentSetRV):
    def __init__(self,currGene,valRange,currNetwork,allNetworks,weightMap,temperature,h):
        super(ParentSetRVFw1,self).__init__(currGene,valRange,currNetwork,allNetworks,weightMap,temperature,h);
        
    def conditionalProbability(self,value, rvList):
        self.setCurrentValue(value);
        wlKernel = WLKernel([self.currNet, self.allNets[-1]]);
        diffScore = wlKernel.computeWLKernelFramework1(self.h); 
        return [diffScore, self.weightMap.getWeight(self.currGene, self.getCurrentValue())];

#Random variable for hyper nodes in 1st framework
class ParentSetStarRV(ParentSetRV):
    def __init__(self, currGene, valRange, currNetwork, allNetworks, temperature,h):
        super(ParentSetStarRV, self).__init__(currGene, valRange, currNetwork, allNetworks, None, temperature,h);
    def conditionalProbability(self,value,rvList):
        self.setCurrentValue(value);
        wlKernel = WLKernel(self.allNets);
        diffScore = wlKernel.computeWLKernelFramework1(self.h); 
        return diffScore;    
    def normalizeDistribution(self,db):
        scores = self.distribution;
        mmx = np.amax(scores);
        myE = np.subtract(mmx,scores);
        temp = np.exp(np.multiply(-1.0*self.beta,myE));
        sums = np.sum(temp);
        self.distribution = [t*1.0/sums for t in temp]; 

def discreteSample(distribution, c, temp):
    counter = 0;
    tot = 0;
    for d in distribution:
        if tot + d > c:
            break;
        else:
            tot = tot + d;
            counter = counter + 1;
    return temp[counter];
        
