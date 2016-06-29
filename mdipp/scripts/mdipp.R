#!/usr/bin/Rscript

library(argparse)
source('/scripts/analysis.R')

#argument parsing - really ugly as the C++ code wants its arguments without dashes
parser = ArgumentParser()
parser$add_argument('argstr', nargs='+')
args = parser$parse_args()

#so, which files do we actually have? load them
flags = c(0,0,0,0)
for (i in 1:(length(args$argstr)-1))
{
	if (args$argstr[i] == 'N')
	{
		yg = loadDataGauss(args$argstr[i+1])
		flags[1] = 1
	}
	if (args$argstr[i] == 'GP')
	{
		ygp = loadDataGP(args$argstr[i+1])
		flags[2] = 1
	}
	if (args$argstr[i] == 'M')
	{
		ymn = loadDataMultinom(args$argstr[i+1])
		flags[3] = 1
	}
	if (args$argstr[i] == 'BW')
	{
		ybw = loadDataBagOfWords(args$argstr[i+1])
		flags[4] = 1
	}
}
#this may well be the ugliest workaround code I have ever written
#as something goes to hell when I append the data structures to a list
if (identical(flags,c(1,1,1,1)))
{
	datafiles = list(yg,ygp,ymn,ybw)
}
if (identical(flags,c(1,1,1,0)))
{
	datafiles = list(yg,ygp,ymn)
}
if (identical(flags,c(1,1,0,1)))
{
	datafiles = list(yg,ygp,ybw)
}
if (identical(flags,c(1,0,1,1)))
{
	datafiles = list(yg,ymn,ybw)
}
if (identical(flags,c(0,1,1,1)))
{
	datafiles = list(ygp,ymn,ybw)
}
if (identical(flags,c(1,1,0,0)))
{
	datafiles = list(yg,ygp)
}
if (identical(flags,c(1,0,1,0)))
{
	datafiles = list(yg,ymn)
}
if (identical(flags,c(0,1,1,0)))
{
	datafiles = list(ygp,ymn)
}
if (identical(flags,c(1,0,0,1)))
{
	datafiles = list(yg,ybw)
}
if (identical(flags,c(0,1,0,1)))
{
	datafiles = list(ygp,ybw)
}
if (identical(flags,c(0,0,1,1)))
{
	datafiles = list(ymn,ybw)
}
if (identical(flags,c(0,0,0,1)))
{
	datafiles = list(ybw)
}
if (identical(flags,c(0,0,1,0)))
{
	datafiles = list(ymn)
}
if (identical(flags,c(0,1,0,0)))
{
	datafiles = list(ygp)
}
if (identical(flags,c(1,0,0,0)))
{
	datafiles = list(yg)
}

#making use of Sam's functions here
mcmc = readMdiMcmcOutput("mcmc_chains.csv")
cpsm = generateConsensusPSM(tail(mcmc))
nclust = apply(getClustersOccupied(tail(mcmc,101)),2,median)
cp = extractPSMClustPartition(cpsm, nclust, datafiles)
write.csv(cp, "cluster_partition.csv")

#make PSM plot
png('PSM_plot.png',width=6,height=6,units='in',res=300)
par(mar=c(1,2,1,1)); plotConsensusPSM(cpsm, datafiles, median(nclust), ann = TRUE)
dev.off()

#dig out the PSM plot consensus clusters and export them
clusdirt = cutree(cpsm$hc,median(nclust))[cpsm$hc$ord]
#ordered as on the PSM plot
clus_nums = rev(unique(clusdirt))
lines = c()
for (i in clus_nums)
{
	genes = names(clusdirt[clusdirt==i])
	header = paste('---Cluster ',toString(i),' ---',sep='')
	#reverse the genes as they're listed from the bottom of the plot
	lines = c(lines,header,rev(genes))
}
fileConn<-file("clusters_PSM_plot.txt")
writeLines(lines,fileConn)
close(fileConn)

#ordered by cluster number
clus_nums = sort(unique(clusdirt))
lines = c()
for (i in clus_nums)
{
	genes = names(clusdirt[clusdirt==i])
	header = paste('---Cluster ',toString(i),' ---',sep='')
	#reverse the genes as they're listed from the bottom of the plot
	lines = c(lines,header,rev(genes))
}
fileConn<-file("clusters.txt")
writeLines(lines,fileConn)
close(fileConn)