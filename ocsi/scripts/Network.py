import numpy as np;
class Vertex(object):
    #here, for directed networks, neighbours are the child nodes
    #for undirected networks, parents is nulls
    #label is an integer that marks the same orthologues hopefully starting from 0
    #neighbours and parents are both lists of vertices
    #parents and neighbours must not be None
    def __init__(self, label, neighbours, parents):
        self.originalLabel = label;
        self.label = label;
        self.temp = '';
        self.neighbours = neighbours;
        self.parents = parents;
    
    #getters and setters for numerical labels     
    def getLabel(self):
        return self.label;
    def setLabel(self, l):
        self.label = l;
        
    #reset label should be called after running kernels that might mess up label
    #it resets label to its original values when the node is initialized
    def resetLabel(self):
        self.label = self.originalLabel;
        
        #getters and setters for parent list, only replacing the list when set
    def getParents(self):
        return self.parents;
    def setParents(self, newParents):
        self.parents = newParents;
    
    #getters and setters for neighbour list, only replacing the list when set
    def getNeighbours(self):
        return self.neighbours;
    def setNeighbours(self, neighbours):
        self.neighbours = neighbours;
    
    #get all neighbours' label, sort them in numerical value, convert to string and concatenate all 'string labels'
    def getNeighbourLabels(self):
        nlabels = [vertex.getLabel() for vertex in self.getParents()] + [self.getLabel()];
        #nlabels = [vertex.getLabel() for vertex in self.getParents()];
        nlabels.sort();
        return ''.join('g'+str(l) for l in nlabels);
        
    
    #make temporary label (for used in WL kernel) by concatenate current label (converted into string)
    #with neighbours sorted 'string label'
    def updateTempLabel(self):
        self.temp = str(self.getLabel()) + self.getNeighbourLabels();
        return self.temp;
    def getTempLabel(self):
        return self.temp;
    def setTempLabel(self, t):
        self.temp = t;
    def setTempToLabel(self):
        self.label = self.temp;
    def __str__(self):
        return str(self.originalLabel);            

class GRNetwork:
    #indlist is incident list, which is only a list of vertices
    #incidency is indicated from within each of the vertex object itself by parents and neighour lists, not by the graph object
    def __init__(self, indlist):
        self.indlist = indlist;
    
    #getter for incident list
    def getIndList(self):
        return self.indlist;
        
    #add a directed edge into the graph
    def addDirectedEdge(self,v1,v2):
        v1.getNeighbours().append(v2);
        v2.getParents().append(v1);
        
        #if (v1,v2) is an edge, remove it by removing v2 from v1's neighbour and removing v1 from v2's parent
    def removeDirectedEdge(self, v1, v2):
        if v2 in v1.getNeighbours() and v1 in v2.getParents():
            v1.getNeighbours().remove(v2);
            v2.getParents().remove(v1);    
    
    #print labels
    def printLabels(self):
        allLabels = [vertex.getLabel() for vertex in self.indlist];
        print(allLabels);
        
    #replace parent set of a node by a new set by changing adding/removing edges
    def updateParent(self, child, newParents):
        for p in child.getParents():
            p.getNeighbours().remove(child);        
        child.setParents([]);
        for p in newParents:
            self.addDirectedEdge(p, child);
    
    #initialize parental set, should only be used once if the graph is to be initialized programmatically
    def initializeParent(self,child,newParents):
        for p in newParents:
            self.addDirectedEdge(p, child);
    
    #reset labels of all vertices in the graph into their original labels (to be called after computing WL kernel)
    def resetLabels(self):
        for v in self.indlist:
            v.resetLabel();

    #output graphs
    def __str__(self):
        rep = '';
        for v in self.indlist:
            rep = rep + str(v)+"\n-->";
            for n in v.getNeighbours():
                rep = rep + str(n)+",";
            rep = rep + "\n<--";
            for n in v.getParents():
                rep = rep + str(n)+",";
            rep = rep + "\n------\n";
        rep = rep + "--------------------------------\n";
        return rep;
    
         
           

        
                
