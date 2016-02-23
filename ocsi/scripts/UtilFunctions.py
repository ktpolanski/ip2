from Network import GRNetwork;
import numpy as np;
from OCSIUtilities import discreteSample as mysample;
from itertools import combinations as combinations;
import pdb;
#build keys of the form #Species#GeneID[#ParentID1#ParentID2...]
def buildKey(species, gid, pas):
    k = str(species)+'#'+gid;
    temp = [];
    if len(pas.strip()) > 0:
        for p in pas.split(':'):
            temp.append('#'+p);
    temp.sort();
    k += ''.join(temp);
    return k;

#build a network from a gene dictionary and parent set distribution (using negative log likelihood)
def buildNetwork(gdict,paDict, weightMap ):
    net = GRNetwork( list(gdict.values()) ); 
    for k,v in gdict.items():
        allPaSets = paDict[k];
        distribution = [];
        for paSet in allPaSets:
            speciesId = v.getSpeciesId();
            gid = v.getGeneId();
            pasId = ':'.join([pa.getGeneId() for pa in paSet]);
            distribution.append(weightMap[buildKey(speciesId, gid, pasId)]);
        initPa = mysample(distribution, np.random.random() , allPaSets);
        net.initializeParent(v, initPa);
    return net;
def buildNetwork2(gdict, paDict, logWeightMap):
    net = GRNetwork( list(gdict.values()) ); 
    for k,v in gdict.items():
        allPaSets = paDict[k];
        dist1 = [];
        for paSet in allPaSets:
            hkey = buildHKey(v,paSet);
            dist1.append(logWeightMap[hkey]);
        dist2 = np.array(dist1);
        distribution = dist2/np.sum(dist2);
        initPa = mysample(distribution, np.random.random() , allPaSets);
        net.initializeParent(v, initPa);
    return net;

def parseNetwork(gdict, netfile):
    net = GRNetwork ( list(gdict.values()) );
    for line in netfile:
        temp = line.strip();
        if temp == '<end>':
            break;
        elif temp == '<begin>' or temp == '':
            continue;
        else:
            parts = temp.strip().split(',');
            currGene = gdict[parts[0].strip()];
            paSet = [gdict[gid.strip()] for gid in parts[1:]];
            net.updateParent(currGene, paSet);
    return net;


def buildHKey(geneChild, genePaList):
    hkey = str(geneChild.getLabel());
    temp = ["#"+str( p.getLabel() ) for p in genePaList];
    temp.sort();
    return hkey+''.join(temp);

class HyperNetworkFactory:
    def __init__(self, allHyperNodes, weightMap, depth):
        self.nodes = allHyperNodes.values();
        self.paStore = [];
        self.depth= depth;
        self.weightMap = weightMap;
        for n in self.nodes:
            self.paStore.append([]);
    def buildPaStore(self):
        for i in range(len(self.nodes)):
            indices = [];
            for d in range(self.depth + 1):
                indices += combinations(range(i)+range(i+1,len(self.nodes)), d);
            self.paStore[i] = [ [self.nodes[ii] for ii in pasetInd] for pasetInd in indices ];
    def buildPaDistribution(self, index):
        dist1 = [];
        for pas in self.paStore[index]:
            hkey = buildHKey(self.nodes[index], pas);
            if (hkey in self.weightMap):
                dist1.append (self.weightMap[hkey]);
            else:
                dist1.append(0);
        dist2 = np.array(dist1);      
        distribution = dist2/np.sum(dist2);
        return distribution;
    def buildNetwork(self):
        net = GRNetwork(self.nodes);
        c = np.random.random(len(self.nodes));
        self.buildPaStore();
        for index in range(len(self.nodes)):
            ntempDist = self.buildPaDistribution(index);
            pachoice = mysample(ntempDist, c[index], self.paStore[index]);
            net.initializeParent(self.nodes[index], pachoice);
        return net;
