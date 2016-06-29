import argparse
import sys
import os

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument(dest='file', type=argparse.FileType('r'))
	args = parser.parse_args()
	
	#white space weirdness for days - .readlines() fails, but .read().splitlines() works
	csv_guts = args.file.read().splitlines()
	args.file.close()
	
	#this thing will exist in the file, detailing how we're doing
	f = open('scores.txt','r')
	probes_guts = f.readlines()
	f.close()
	
	probes = []
	for line in probes_guts:
		blaa = line.split()
		probes.append(blaa[0])
	
	f = open(args.file.name,'w')
	#double header days
	f.write(csv_guts[0])
	f.write('\n'+csv_guts[1])
	csv_guts2 = csv_guts[2:(len(csv_guts)+1)]
	#toggle to assess if we're making GP2S progress or just looping infinitely because of something else
	toggle=0
	for line in csv_guts2:
		blaa = line.split(',')
		if blaa[0] not in probes:
			f.write('\n'+line)
		else:
			#we've got a thing for which we got a GP2S score since last time. we're good
			toggle=1
	f.close()
	
	#if the toggle is still at 0, nothing got computed since the previous check
	#and that may imply something fishy is going on. let's make it fail three times though
	if toggle==0:
		if os.path.isfile('dummy.txt'):
			f = open('dummy.txt','r')
			line = f.readline()
			f.close()
			line = int(line)
		else:
			#no dummy file created, so this is the first time the script is running
			#and, as such, automatically the first crash it encountered
			line = 0
		#if the counter is on 2, then this is the third failure
		if line==2:
			os.remove('dummy.txt')
			sys.exit(1)
		#otherwise, just let this know it failed before
		f = open('dummy.txt','w')
		f.write(str(line+1))
		f.close()
	else:
		#all clear!
		f = open('dummy.txt','w')
		f.write('0')
		f.close()

if __name__== '__main__':
	main()