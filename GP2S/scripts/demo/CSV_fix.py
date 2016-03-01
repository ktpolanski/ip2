import argparse

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
	for line in csv_guts2:
		blaa = line.split(',')
		if blaa[0] not in probes:
			f.write('\n'+line)
	f.close()

if __name__== '__main__':
	main()