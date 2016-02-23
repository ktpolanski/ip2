from numpy import dot;     
import pdb;
class WLKernel:
    def __init__(self, graphList):
        self.graphs = graphList;
        self.featureWLTrees = [[] for g in self.graphs];
        self.counter = max([max([v.getLabel() for v in g.getIndList()]) for g in self.graphs]) + 1;
        self.newAlphabet = [];
        
    #using sets, slow! update to be described algorithm in the paper for speed  latter on
    def wlUpdate(self):
        newLabelMap = {};
        self.newAlphabet = [];
        for g in self.graphs:
            for vertex in g.getIndList():
                temp = vertex.updateTempLabel();
                if not(temp in newLabelMap):
                    self.newAlphabet.append(self.counter);
                    newLabelMap[temp] = self.counter;
                    self.counter = self.counter + 1;
                vertex.setTempLabel(newLabelMap[vertex.getTempLabel()]);
        for g in self.graphs:
            for v in g.getIndList():
                v.setTempToLabel();
            
    def countfrequency(self, graph, alphabet):
        allLabels = [vertex.getLabel() for vertex in graph.getIndList()];
        freq = [allLabels.count(letter) for letter in alphabet];
        return freq;
        
    def wlRecordFeature(self):
        alphabet = self.newAlphabet;            
        #alphabet.sort();
        i = 0;
        while i < len(self.graphs):
            self.featureWLTrees[i] = self.featureWLTrees[i] + self.countfrequency(self.graphs[i],alphabet); 
            i = i + 1;   
        #pdb.set_trace();
        
    def wlInitFeature(self):
        self.newAlphabet = list(range(self.counter));
        self.wlRecordFeature();
    
    #compute the WL kernel by dot product of the feature trees and reset all labels of the graphs
    def computeWLKernelFramework2(self, h):
        self.wlInitFeature();
        for i in range(h):
            self.wlUpdate();            
            self.wlRecordFeature();
        for g in self.graphs:
            g.resetLabels();
        score =  0;
        i = 0;
        while i < len(self.graphs):
            j = i + 1;
            while j < len(self.graphs):
                score = score + dot(self.featureWLTrees[i], self.featureWLTrees[j]);
                j = j + 1;
            i = i + 1; 
        return score;
    def computeWLKernelFramework1(self, h):
        self.wlInitFeature();
        for i in range(h):
            self.wlUpdate();            
            self.wlRecordFeature();
        for g in self.graphs:
            g.resetLabels();
        score =  0;
        i = 0;
        while i < len(self.graphs) - 1:
            score = score + dot(self.featureWLTrees[i], self.featureWLTrees[-1]);
            i = i + 1; 
        return score;

