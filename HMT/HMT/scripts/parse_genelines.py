import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('findstr',type=str)
args = parser.parse_args()

#equal sign support!
findstr = args.findstr
if findstr[-1]!='=':
	findstr+='='

#import using fancy pandas technology
#(and then immediately get rid of the fancy pandas technology)
genelines = pd.read_csv('genelines.gff3', sep='\t', index_col=None, header=None)
genelines = genelines.values

#attempt to pull out the appropriate gene name and clump it into the final field
#(instead of the whole ;-separated shebang)
error_toggle = 0
for i in range(genelines.shape[0]):
	line = genelines[i,-1].split(';')
	local_toggle = 0
	for field in line:
		if findstr == field[:len(findstr)]:
			genelines[i,-1] = field[len(findstr):]
			local_toggle = 1
			break
	if local_toggle == 0:
		error_toggle = 1
		sys.stderr.write('Failed to find specified header ('+findstr[:-1]+') in GFF3 line: '+'\t'.join(genelines[i,:]))

#if our input is borked, we can't process it
if error_toggle == 1:
	sys.exit(1)

#okay then. have the processed thing ready. need to format it accordingly
genelines = genelines[:,[0,3,4,8,3,6]]
#for reasons I do not understand, the next to last column of the BED has to be a number.
for i in range(genelines.shape[0]):
	genelines[i,-2] = 1

#that's it. save it
np.savetxt('genelines.bed',genelines,fmt='%s\t%i\t%i\t%s\t%i\t%s')