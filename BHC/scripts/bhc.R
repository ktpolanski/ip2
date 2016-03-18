#!/usr/bin/Rscript

suppressPackageStartupMessages(require(argparse))
suppressPackageStartupMessages(require(BHC))

#argument parsing
parser = ArgumentParser()
parser$add_argument('file', nargs=1, help='CSV with gene expression, first column gene names, first row time points')
parser$add_argument('--Mode', dest='mode', default='mutinomial', help='Mode of operation (multinomial/time-course). Default: multinomial')
args = parser$parse_args()

#basic data prep
data = read.csv(args$file,header=TRUE,row.names=1,check.names=FALSE)
genes = rownames(data)
samples = colnames(data)
data = data.matrix(data)

if (args$mode == 'time-course')
{
	#get them time points as time points
	samples = as.integer(samples)
}