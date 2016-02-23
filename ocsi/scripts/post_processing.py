import pdb;
import os, glob, shutil;
def process_line(line,burnin):
    parts = line.split(';');
    myArray = parts[burnin:];
    adict = {};
    for paset in myArray:
        if paset in adict:
            adict[paset] += 1;
        else:
            adict[paset] = 1;
    mysum = 0;
    for k,v in adict.items():
        mysum += v;
    for k,v in adict.items():
        adict[k] = 1.0*v/mysum;
    return adict;

def buildMatrix(filename, burnin):
    matrix = {};
    rfn = filename+"RawMCMC";
    with open(rfn,'r') as f:
        for line in f:
            parts = line.strip().split(':');
            speciesId = int(parts[0]);
            geneId = parts[1];
            if not(speciesId in matrix):
                matrix[speciesId] = {};
            matrix[speciesId][geneId] = process_line(parts[2],burnin);
    return matrix;

def writeOutput(filename,matrix):
    firstTime = True;
    for k,v in matrix.items():
        speciesname = str(k);
        if speciesname == '-1':
            speciesname = 'hypernet';
        outputName = filename+'OcsiWeight'+speciesname+".csv";
        fw = open(outputName,'w');
        for k1 in sorted(matrix[k].keys()):
            for k2,v2 in matrix[k][k1].items():
                paset = k2.replace(',',':');
                fw.write(k1+','+paset+','+repr(v2));
                fw.write('\n');
        fw.close();
        makeMarginalMatrix(outputName,filename+'Marginal'+speciesname+'.csv');
        #move the OCSI weight file into a utility folder
        utilFolder = filename+"UtilFolder";
        if not os.path.exists(utilFolder):
            os.makedirs(utilFolder);
        elif firstTime:
            firstTime = False;
            filelist = glob.glob(os.path.join(utilFolder,"*"));
            for f in filelist:
                os.remove(f);
        shutil.move(outputName,utilFolder);

def makeMarginalMatrix(fin,fout):
    mydict = {};
    with open(fin,'r') as f:
        for line in f:
            parts = line.strip().split(',');
            curGene = parts[0];
            pasets = parts[1];
            palist = pasets.split(':') 
            weight = float(parts[2]);
            if not(curGene in mydict):
                mydict[curGene] = {};
            for p in palist:
                if p in mydict[curGene]:
                    mydict[curGene][p] += weight;
                else:
                    mydict[curGene][p] = weight;
    idlist= [];
    try: #this is just to work on dream data, usually this will crash, but does not matter if it works
        idlist = sorted(mydict.keys(),key=lambda x: int(x[1:]));
    except: #if the above crashes, then just do this
        idlist = mydict.keys();
    fw = open(fout,'w');
    fw.write(','+','.join(idlist)+'\n');
    for gid in idlist:
        temp = [gid];
        for gid2 in idlist:
            if gid2 == gid or not(gid2 in mydict[gid]):
                temp.append('0');
            else:
                temp.append(repr(mydict[gid][gid2]));
        fw.write(','.join(temp)+'\n');
    fw.close();
            
            
            
        
