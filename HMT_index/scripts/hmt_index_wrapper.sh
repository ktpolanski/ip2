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
# $8 - toggle whether to remove promoter overlapping bits with gene sequences
# $9 - toggle whether to do Uniprobe to MEME conversion
# $10 - toggle whether to include 5' UTR sequence

#start off by filtering the .gff3 to gene lines only
if [ ! -f annot.gff3 ]
then
	cp $2 annot.gff3
fi
grep -P '\tgene\t' annot.gff3 > genelines.gff3

#strip the potential FASTA line breaks. creates genome_stripped.fa
if [ ! -f genome.fa ]
then
	cp $1 genome.fa
fi
python3 /scripts/strip_newlines.py

#create the .genome file
samtools faidx genome_stripped.fa
cut -f 1-2 genome_stripped.fa.fai > bedgenome.genome

#parse up the .bed for promoter extraction
python3 /scripts/parse_genelines.py $3
#the python script takes the genelines.gff3 file and makes a genelines.bed out of it
bedtools flank -l $4 -r 0 -s -i genelines.bed -g bedgenome.genome > promoters.bed
#remove overlapping promoter chunks
if [ $8 == '--NoOverlap' ]
then
	bedtools subtract -a promoters.bed -b genelines.bed > promoters2.bed
	mv promoters2.bed promoters.bed
fi
#possibly add 5' UTR
if [ $10 == '--UseUTR' ]
then
	python3 /scripts/parse_utrs.py
fi
python3 /scripts/parse_promoter_lengths.py
bedtools getfasta -fi genome_stripped.fa -bed promoters.bed -s -fo promoters_rough.fa
#this results in some really crappy nomenclature for gene names
#so let's make promoters.fa ourselves
python3 /scripts/parse_promoters.py

#now we can actually FIMO our way to victory
/root/meme/bin/fasta-get-markov promoters.fa > promoters.bg
#FIMO barfs ALL the output. that's not good. time for individual FIMOs
#on individual MEME-friendly motif files too
mkdir memefiles
mkdir logos
#optionally turn Uniprobe into MEME. MEME sucks as a format, Uniprobe is easier
if [ $9 == '--Uniprobe' ]
then
	/root/meme/bin/uniprobe2meme $5 > MEME-motifs.txt
	python3 /scripts/parse_memefile.py MEME-motifs.txt
else
	python3 /scripts/parse_memefile.py $5
fi
#now loop over all the individual motif files
for fid in memefiles/*.txt
do
	bfid=$(basename ${fid/.txt/})
	/root/meme/bin/fimo --text --thresh $6 --verbosity 1 --bgfile promoters.bg $fid promoters.fa > fimo.txt
	#this appends to fimo_found.txt
	python3 /scripts/parse_matrix.py $7 $6
	#generate logo too while we're at it
	/root/meme/bin/ceqlogo -i $fid -m 1 -o logos/$bfid.eps
	#convert to PNG because Paul says so, via ImageMagick because reasons
	convert logos/$bfid.eps logos/$bfid.png
	rm logos/$bfid.eps
done

#there's a lot of intermediate files that need blanking
rm -r memefiles
if [ $9 == '--Uniprobe' ]
then
	rm MEME-motifs.txt
fi
rm bedgenome.genome
rm genelines.bed
rm genelines.gff3
rm genome.fa
rm genome_stripped.fa
rm genome_stripped.fa.fai
rm promoters.bed
rm promoters.bg
rm promoters.fa
rm promoters_rough.fa
rm promoter_lengths.txt