# Causal Structural Inference (CSI) #

## Introduction ##

CSI is a tool for inferring causal relationships between genes, given
timeseries of expression profiles.  The method was initially developed
and applied to the analysis of expression data in:

    Penfold, Christopher A., and David L. Wild. "How to infer gene
    networks from expression profiles, revisited." Interface Focus 1,
    no. 6 (2011): 857-870.  doi:10.1098/rsfs.2011.0053

The method attempts to discover which genes best explain the changes
in expression of other genes.  Slightly more formally, CSI infers the
probability that any given *target* gene is being regulated by any
*parent* gene.  These probabilities can then be employed to create a
network graph describing the inferred regulatory network.

Instead of working directly with parents, CSI instead evaluates the
probabilities of a large number of *parental sets* (i.e. each parental
set contains some number of  parents) predicting the target gene.  The
probabilities associated with each parental sets can then be summed to
calculate the probabilities of any given gene having a causal effect
on another gene.

# File Formats #

The input to CSI is a single CSV file that contains expression data.
The data should have two header rows, the first describing the
replicate/perturbation of this timeseries and the second row
describing the times at which the data was sampled.  The remainder of
the file is the *standardised* expression data, by this it is meant
that the mean of the data is approximately zero and standard deviation
is approximately one.  Depending on your data it may be better to
ensure that each gene has been independently standardised and other
analyses may benefit from standardisation over all genes.

CSI outputs a HDF5 file whose structure is designed to capture all
intermediate information from the analysis.  Other tools able to
translate this into more immediately useful files, such as a HTML
based analysis tool for of the results visualisation.

### Definitions ###

CSV
 ~ Comma separated variables, most commonly available tools will be
   able to create an appropriately formatted CSV data file

HDF5
 ~ Hierarchical Data Format v5, is a binary format output by default
   by CSI.  The format can be directly read by a number of common
   programming languages, including Matlab, R and Python.

## Limitations ##

CSI evaluates every possible parental sets


# Glossary #

Target
 ~ A single gene or other experimentally measured 

Parental Set
 ~ Definition 2a

Power Set
 ~ a

Truncated Power Set
 ~ b
