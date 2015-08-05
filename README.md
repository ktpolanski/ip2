# iPlant Integration Project from Warwick's System Biology #

This repository contains the tools that will be ported into the iPlant
project.  As a starting point, the tools consist of:

* GP2S: A robust Bayesian two-sample test for detecting intervals of
  differential gene expression in microarray time series.
  doi:10.1089/cmb.2009.0175

* CSI: Causal structure identification algorithm for inferring gene
  regulatory networks.  How to infer gene networks from expression
  profiles, revisited. doi:10.1098/rsfs.2011.0053

* HCSI: Hierarchical variant of the CSI.  Nonparametric Bayesian
  inference for perturbed and orthologous gene regulatory
  networks. doi:10.1093/bioinformatics/bts222

* Wigwams: identifying gene modules co-regulated across multiple
  biological conditions. doi:10.1093/bioinformatics/btt728

* Apples: Analysis of Plant Promoter-Linked Elements.  http://www2.warwick.ac.uk/fac/sci/systemsbiology/staff/ott/tools_and_software/apples/

* Gradient Tool: not sure where this is from.  Looks for significant
changes in gradient of timeseries.  http://www.plantcell.org/content/23/3/873.long


* VBSSM: A Bayesian approach to reconstructing genetic regulatory
  networks with hidden factors.  doi:10.1093/bioinformatics/bti014

* MVBSSM: not sure how this variant alters the above.

* MemeLab: motif analysis in clusters.
  doi:10.1093/bioinformatics/btt248

* Wellington: a novel method for the accurate identification of
  digital genomic footprints from DNase-seq data.
  doi:10.1093/nar/gkt850

Some of these projects are already written in languages compatible
with the iPlant infrastructure, but others are based in Matlab and
will be ported to Python and C++ for compatibility and performance.

# TODO #

# General Notes #

Tools could have a header at the beginning of their main report/PDF
containing the name of the creating tool along with citation
information.

## iPlant / Agave ##

The file system used within iPlant is called *iRODS*.  It comes with
its own set of *iCommands* for performing `ftp` like operations at the
command line.  Alternatively there are GUI apps like iDrop and a
`fuse` plugin available.  There are also nice Python libraries
available for talking to iRODS servers, for example [python-irodsclient].

To connect to the iPlant servers, you need the following details:

    host: data.iplantcollaborative.org
    port: 1247
    zone: iplant
    username: <iplant username>
    password: <iplant password>

## Python Dependencies ##

These are the minimal Python libraries that would need to be installed:

    pip install numpy scipy pandas matplotlib GPy

In order to run an interactive exploration in iPython, these are also
needed [google]:

    pip install ipython pyzmq jinja2 tornado jsonschema

[python-irodsclient]: https://github.com/iPlantCollaborativeOpenSource/python-irodsclient
