#!/bin/bash
set -e

#inputs:
# $1 - .gff3 annotation, as supplied to the indexing function
# $2 onwards - python hypergeometric test arguments

grep -P '\tgene\t' $1 > genelines.gff3