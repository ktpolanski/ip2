#!/bin/bash
set -e

#inputs:
# $1 - .gff3 annotation, as supplied to the indexing function
# $2 - .gff3 header to use, as supplied to the indexing function (TAIR 'gene_id=')
# $3 onwards - python hypergeometric test arguments

#start off similarly to the indexing wrapper - get out the gene lines
grep -P '\tgene\t' $1 > genelines.gff3
python3 /scripts/parse_genelines.py $2
#now that we have genelines.bed, we can pull out the gene IDs to use as our universe
cut -f 4 genelines.bed > universe.txt

#and that's all, folks. HMT time
python3 /scripts/hmt.py --Universe universe.txt ${@:3}

#placeholder file cleanup to make iPlant output cleaner
rm genelines.bed
rm genelines.gff3
rm universe.txt