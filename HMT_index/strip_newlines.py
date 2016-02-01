fid = open('genome.fa','r')
fid2 = open('genome_stripped.fa','w')
toggle = ''
for line in iter(fid):
	#lose the newline, whatever its source
	line = line.rstrip()
	#sometimes there's just newline
	if not line:
		continue
	if line[0]=='>':
		#we got a gene ID on our hands
		line = toggle+line+'\n'
		toggle = '\n'
	print(line,end='',file=fid2)
fid.close()
fid2.close()