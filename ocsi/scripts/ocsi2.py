import argparse;
from AbstractGibbs import AbstractGibbs as AbstractGibbs;
from OCSIUtilities import ParentWeightMap as ParentWeightMap;
from OCSIUtilities import ParentSetRVFw1 as ParentSetRVFw1;
from OCSIUtilities import ParentSetRVFw2 as ParentSetRVFw2;
from OCSIUtilities import ParentSetStarRV as ParentSetStarRV;
from OCSIUtilities import Gene as Gene;
from UtilFunctions import buildKey as buildKey;
from UtilFunctions import buildNetwork as buildNetwork;
from UtilFunctions import buildNetwork2 as buildNetwork2;
from UtilFunctions import parseNetwork as parseNetwork;
from UtilFunctions import buildHKey as buildHKey;
from UtilFunctions import HyperNetworkFactory as HyperNetworkFactory;
import math;
import pdb;
import pandas as pd;
from post_processing import buildMatrix as buildOutputMatrix;
from post_processing import writeOutput as writeOutput;
################ Dummy class & utilities functions for parsing######
class C:
    pass;
################ Argument parsing ##################################
parser = argparse.ArgumentParser();
parser.add_argument('-i', required = True, nargs = '+', \
                        help = 'CSI output files to be used by OCSI');
parser.add_argument('-beta', required = True,  type = float, help = 'Temperatures');
parser.add_argument('-fw', default = 1,  type = int, help = 'Ocsi framework 1 or 2');
parser.add_argument('-l',  required = True, \
                        help = 'Labels for genes, in the order that they apper in CSI output files');
parser.add_argument('-br', default = 100, type = int, help = 'Burnin');
parser.add_argument('-d', default = 2, type = int, help = 'Size of parental set of hyper network (must be consistent with CSI)');
parser.add_argument('-fn', help = 'A file containing network topology (incident list) of species with known network topology');
parser.add_argument('-g',default = 1000, type = int, help = 'Number of steps for Gibbs sampling');
parser.add_argument('-o', default = 'out', help = 'Output files\' name prefix for OCSI');
parser.add_argument('-wlh', default = 1, type = int, help = 'WL kernel depth');
parser.add_argument('-init', default = 0, type = int, help = 'Selecting initialization methods for networks');
myargs = C();
parser.parse_args(namespace = myargs);

################### Reading parsed arguments ################################
beta = myargs.beta;  #temperature
framework = myargs.fw; #framework 1 uses hypernetwork and framework 2 use direct pairwise network leveraging
csiFileNames = myargs.i; #input files from csi-output
labelFile = myargs.l; #files that contains labels
fixedNetFileName = myargs.fn; #files that contain network structures of known species
labels = []; #empty array to store all labels of all genes of all species
depth = myargs.d; #parental depth
gibbsSteps = myargs.g; #how many Gibbs iterations to do
outputFileName = myargs.o; #output file name

h = myargs.wlh; #how deep WL go into
burnin = myargs.br;# how many burnin step
myinit = myargs.init;
if myinit == 0:
    myinit = framework;


#################################### an exhaustive list of important variables used in the driver ########################
allHyperNodes = {}; #if using framework 1: dictionary of hyper nodes. keys are hyper node labels. labels == gene id for hyper nodes
logweightSum = {}; #MAJOR MISLEADING NAME, IT'S NOT LOG. "summary weight" i.e. the product weight of corresponding genes (same labels) in each species
allGenes = []; #list of species, each item of the list is a dictionary of all genes based on geneID
listAllGenes = [];#list of species, each item of the list is a list of all genes in the order they appear in csi output, essentially
# this is the same as allGenes, but elements are not hashmap but list
csiWeight = {}; #weightmap that takes species,target gene and parent to map into csi weight
shouldUpdate = []; #list indicate if a species network is known or not
paStore = []; # a list of dictionary where the ith item in the list corresponds to the ith species. each key in the dictionary of the ith species is the id of a gene in that species
allNets = []; # a list of all GRN if framework 1 is used, then the hyper net is always the last network in this list
rvs = [];# list of random variable (parental sets object for genes) to be used in Gibbs sampler
idToLabelMap = {}; #a hashmap that takes gene id to map into orthologue label
noDistinctLabel = 0;
######################## read in label file ###################
lbr = open(labelFile,'r');
lineno = 1;
for line in lbr:
    gidList = line.strip().split(',');
    if len(gidList) > 0:
        noDistinctLabel += 1;
        for gid in gidList:
            t = gid.strip();
            idToLabelMap[t] = lineno;
    lineno += 1;
lbr.close();
############################################ read data from csi ###################
speciesId = 0;
for fname in csiFileNames:
    print("Parsing input from "+fname);
    gdict = {}; paStore.append({}); listAllGenes.append([]);
    df = pd.read_csv(fname,header = None); df.fillna('#', inplace = True);
    #create all genes for this species
    allGeneIds = df[0].unique();
    for gid in allGeneIds:
        gdict[gid] = Gene(idToLabelMap[gid.strip()], gid, speciesId);
        #initialize paStore for this gene
        paStore[speciesId][gid] = [];
        #collect this gene and put it into a list (later used to make random variable object)
        listAllGenes[speciesId] = [gdict[gid]] + listAllGenes[speciesId];
    #read through data frame row by row (i know, it's not the most efficient method)
    for rindx in df.index:
        temp = df.iloc[rindx];
        #first gene id is the target gene object
        target = gdict[temp[0].strip()];
        #list of parents gene object
        parents = [];
        if not (temp[1] == "#"):
            paidlist = temp[1].strip().split(':');
            parents = [gdict[paid.strip()] for paid in paidlist];
        #add parents into parent store
        paStore[speciesId][target.getGeneId()].append(parents);
        #build key from species, target gene and its parents add corresponding csi weight into dictionary
        keys = buildKey(speciesId, target.getGeneId(), ':'.join([p.getGeneId() for p in parents]));
        csiWeight[keys] = float(temp[2]);
        #hyper weight to either initialize hypernet or to initialize using 2nd initialization method
        hkey = buildHKey(target, parents);
        if hkey in logweightSum:
            logweightSum[hkey] *= csiWeight[keys];
        else:
            logweightSum[hkey] = csiWeight[keys];                    
    allGenes.append(gdict);
    speciesId += 1;


####build table lookup of parent weight for later use#######################
csiWeightMap = ParentWeightMap(csiWeight);
#### INITIALIZE network structures################################################
#sample network from csi output
species = 0;
for gdict in allGenes:
    if myinit == 1:
        allNets.append(buildNetwork(gdict, paStore[species], csiWeight));
    elif myinit == 2:
        allNets.append(buildNetwork2(gdict,paStore[species],logweightSum));
    species += 1;

#parse known network topology
if fixedNetFileName is not None:
    networkFile = open(fixedNetFileName,'r');
    geneParser = open(fixedNetFileName,'r');
    tempGDict = {};
    for line in geneParser:        
        if line.strip() == '<end>':
            allNets.append(parseNetwork(tempGDict,networkFile));
            tempGDict = {};
            species += 1;
        else:
            parts = line.strip().split(',');
            targetId = parts[0].strip();
            tempGDict[targetId] = Gene(idToLabelMap[targetId], targetId, species);
    geneParser.close();
    networkFile.close();
#build random variable objects to be passed into gibbs sampling
#ordering the list of RV as follow:
# 1st gene in 1st species, 1st gene in 2nd species, 1st in 3rd species...
# 2nd gene in 1st species, 2nd gene in 2nd species ....
# ....
# nth gene in 1st species, nth gene in 2nd species ...
keepPoping = True;
while keepPoping == True:
    speciesId = 0;
    keepPoping = False;
    for gl in listAllGenes: #loop through all gene lists of all species
        #only create random variable for species that network topology is not known (shouldUpdate[speciesID] = true)
        #if this species gene LIST still have genes, pop out that gene and make a paRV
        if len(gl) > 0:
            g = gl.pop();
            if framework == 1:
                rvs.append(ParentSetRVFw1(g, paStore[speciesId][g.getGeneId()],\
                                       allNets[speciesId], allNets, csiWeightMap, beta, h));
            else:
                rvs.append(ParentSetRVFw2(g, paStore[speciesId][g.getGeneId()],\
                                       allNets[speciesId], allNets, csiWeightMap, beta, h));
        if len(gl) > 0: #if this species still have genes, we need to go at least another round
            keepPoping = True;
        speciesId += 1;


######### if working using framework 1, initialize hypernet, make hyper node parent RV #########
if framework == 1:   #initialize hypernework stuff (paset for each hyper node and net topology)
    uniqueLabels = range(noDistinctLabel);
    for al in uniqueLabels:
        allHyperNodes[(al+1)] = Gene(al+1,al+1,-1);
    hNetFactory = HyperNetworkFactory(allHyperNodes, logweightSum, depth);
    allNets.append(hNetFactory.buildNetwork()); 
    hnodeList = hNetFactory.nodes;
    for gindex in range(len(hnodeList)): #for all hyper nodes in the hyper net, build parental set RV
        temp = hNetFactory.paStore[gindex]; #parent set for hyper node
        #temp.pop(0); temp.append([]); #this is to makesure the empty pa set is at the end (to conform w/ MATLAB)
        rvs.append(ParentSetStarRV(hnodeList[gindex], temp , allNets[-1], allNets, beta, h));

gibbs = AbstractGibbs(rvs);
print("Starting Gibbs sampler:");
gibbs.sample(gibbsSteps);
#####  Saving data ##############
print('Start saving');

sampledValues = gibbs.sampledValues;
fw = open(outputFileName+'RawMCMC','w');
index = 0;
for sv in sampledValues:
    string = [];
    for paSet in sv:
        tstring = '';
        tstring = ','.join([str(n.getGeneId()) for n in paSet]);
        string.append( tstring);
    fw.write(str(rvs[index].currGene.getSpeciesId())+':'+str(rvs[index].currGene.getGeneId())+':');
    fw.write(';'.join(string));
    fw.write('\n');
    index += 1;
fw.close();

outputMatrix = buildOutputMatrix(outputFileName,burnin); #get summary from MCMC chain file
writeOutput(outputFileName,outputMatrix); #create all different output files for

import shutil;
#move raw MCMC file into folder
utilFolder = outputFileName+"UtilFolder";
shutil.move((outputFileName+"RawMCMC"),utilFolder);
fw2 = open(outputFileName+"species.txt","w");
spIndex = 0;
for speName in csiFileNames:
    fw2.write(speName+"    "+str(spIndex)+"\n");
    spIndex += 1;
fw2.close();
shutil.move((outputFileName+"species.txt"),utilFolder);    





            
        
















            
                
