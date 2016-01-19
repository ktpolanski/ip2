#!/bin/bash
set -e

#inputs:
# $1 - genome
# $2 - gff3
# $3 - gff3 header to use (TAIR 'gene_id=')
# $4 - promoter length
# $5 - MEME-friendly motif input
# $6 - significance threshold
# $7 - number of top motifs to take
# $8 - toggle whether to do Uniprobe to MEME conversion
# $9 - toggle whether to remove promoter overlapping bits with gene sequences

#start off by filtering the .gff3 to gene lines only
grep -P '\tgene\t' $2 > genelines.gff3

#create the .genome file
samtools faidx $1
cut -f 1-2 $1.fai > bedgenome.genome

#parse up the .bed for promoter extraction
python3 /scripts/parse_genelines.py $3
#the python script takes the genelines.gff3 file and makes a genelines.bed out of it
bedtools flank -l $4 -r 0 -s -i genelines.bed -g bedgenome.genome > promoters.bed
#remove overlapping promoter chunks
if [ $9 == '--NoOverlap' ]
	then
		bedtools subtract -a promoters.bed -b genelines.bed > promoters2.bed
		mv promoters2.bed promoters.bed
fi
bedtools getfasta -fi $1 -bed promoters.bed -s -fo promoters_rough.fa
#this results in some really crappy nomenclature for gene names
#so let's make promoters.fa ourselves
python3 /scripts/parse_promoters.py

#now we can actually FIMO our way to victory
/root/meme/bin/fasta-get-markov promoters.fa > promoters.bg
#FIMO barfs ALL the output. that's not good. time for individual FIMOs
#on individual MEME-friendly motif files too
mkdir memefiles
#optionally turn Uniprobe into MEME. MEME sucks as a format, Uniprobe is easier
if [ $8 == '--Uniprobe' ]
	then
		/root/meme/bin/uniprobe2meme $5 > MEME-motifs.txt
		python3 /scripts/parse_memefile.py MEME-motifs.txt
	else
		python3 /scripts/parse_memefile.py $5
fi
#now loop over all the individual motif files
for fid in memefiles/*.txt
do
	/root/meme/bin/fimo --text --thresh $6 --verbosity 1 --bgfile promoters.bg $fid promoters.fa >> fimo.txt
	#this appends to fimo_found.txt
	python3 /scripts/parse_matrix.py $7 $4 $6
	#need to wipe fimo.txt to avoid silliness
	rm fimo.txt
done