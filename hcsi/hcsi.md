# hCSI

## The Purpose of the Algorithm

hCSI (hierarchical Causal Structure Inference) is an algorithm designed for inferring regulatory network models from time course datasets showing the same organism responding to different perturbations. The underlying assumption is that the response across the conditions may show some degree of similarity, with elements of a common, shared regulatory network (referred to as a hypernetwork) in use during the individual responses. The distinctness of the stimuli is tuned during the hCSI run, with data for conditions deemed dissimilar having less impact on the hypernetwork and, by extension, other conditions. The algorithm produces network models for each stimulus, as well as an estimate of the hypernetwork.

## How Does It Work?

The basis for hCSI is CSI (Causal Structure Inference), a network inference algorithm designed to propose a regulatory network model for a single time course dataset. CSI operates unter the premise that when a particular parent gene regulates a target (child) gene, then the child's expression profile can be explained by the parent's expression profile at the preceding time point (giving the parent gene time to get its protein to the child gene's promoter and affect its expression). This is formally done by fitting a Gaussian process to a visualisation of the data where one axis is made up of the child's expression and the other axes are time shifted parents. CSI infers parents on a per-child basis, with all eligible parent combinations (up to a specified maximum number of parents) evaluated, and probabilities of each parental set being "the correct one" are computed based on how good that individual fit is. The fits are optimised through Expectation Maximisation, with the values of hyperparameters governing the behaviour of the Gaussian processes tuned to maximise the overall fit quality. CSI was applied to synthetic 10 and 100 gene networks in a comparison with several other algorithms, and performed very well [(Penfold et al., 2011)][penfold2011].

hCSI takes CSI as a base and creates an individual CSI fit for each perturbation's worth of data fed into it. CSI would also be capable of analysing the data, but would produce a single network jointly inferred for all the data at once, whilst hCSI allows for a degree of deviation between the individual stimulus networks. This is achieved through the introduction of the hypernetwork, which can be interpreted as the underlying "master" network that the individual responses are based on, with an additional temperature parameter controlling the level of the hypernetwork's influence on a particular stimulus response. In contrast to CSI and its use of Expectation Maximisation, hCSI makes use of sampling, with a single step involving Gibbs sampling of individual dataset networks and the overarching hypernetwork, and Metropolis-Hastings sampling of the individual dataset hyperparameters and temperatures. Similarly to CSI, hCSI performs these procedures on a per-gene basis, with every gene taking on the role of the child whose expression is to be explained. The underlying mathematics can be found in [Penfold et al., 2012][penfold2012].

## Skippable Technical Minutiae

This section describes the exact details of the implementation, going into more technical detail than the algorithm outline in [Penfold et al., 2012][penfold2012]. It is not necessary for using the iPlant app, but it is here to document the implementation in greater depth.

The model is initialised completely randomly. Each perturbation has a random parental set combination assigned to it as the starting value with equal weighting, and the same happens for the hypernetwork. The starting Gaussian process hyperparameters are sampled from U(0,1) for each dataset independently (in contrast to the original Matlab implementation, where a single set of sampled values are copied over for each dataset as the starting hyperparameters), all of the temperatures are initialised at 0.1.

In contrast to the original Matlab implementation, all four sampling procedures (individual networks, hypernetwork, hyperparameters, temperatures) are performed in tandem, with a default total of 25,000 steps. In Matlab, one of those operations was chosen randomly as a single step, with 100,000 default total steps.

The hyperparameter and temperature sampling steps feature the addition of a random Gaussian variable, which is sampled from N(0,1) and multiplied by a scaling constant. The scaling constant is initialised at 0.1 and re-evaluated every 100 steps, aiming to keep the Metropolis-Hastings acceptance rate at around 0.25. In the case of less than 15 accepted jumps (out of 100), the transition operator is deemed to be not localised enough and the scaling constant is multiplied by 0.9. If over 35 jumps get accepted, the scaling constant is multiplied by 1.1 to make the sampling be less localised. This, as well as the actual sampling of the hyperparameter/temperature values, is performed independently for hyperparameters and temperatures and each dataset.

## Test Run

If you want to take hCSI out for a spin without using your own data, this can be done with the aid of one of the 10-gene synthetic networks originally used in the CSI and hCSI publications. The dataset to be used on input can be found at `ktpolanski/hcsi_testdata/dream4_5.csv` under Community Data. Leave all the parameter values as defaults, except for the process count, which you should set to 10. hCSI is computationally intensive and it may take upwards of an hour for this analysis to be complete.

## Inputs

### Gene Expression CSV

**The only obligatory file you have to provide.** Comma-delimited file, with expression data ordered to have genes as rows and time points as columns. In terms of headers, the first column should contain gene IDs, the first row should contain condition names (repeated for each time point part of the condition), and the second row should contain the corresponding time of the time point in that condition. For reference on formatting, consult `ktpolanski/hcsi_testdata/dream4_5.csv` under Community Data.

**NOTE: hCSI is extremely computationally intensive.** It is greatly recommended to perform a high degree of preliminary analysis and select the most relevant subset of candidate genes to perform hCSI on. Going above 30 genes is really not recommended.

### Parental Set Depth

**Default:** 2

When evaluating parental set combinations, a limitation is put on up to how many parents to sample from the parent pool to create the combinations. As this depth is increased, the number of parental sets to evaluate drastically goes up, making the problem less computationally tractable. ***Increasing this value from 2 is not recommended**, especially if the dataset is larger than the test data provided.

[penfold2011]: http://rsfs.royalsocietypublishing.org/content/1/6/857.short
[penfold2012]: http://bioinformatics.oxfordjournals.org/content/28/12/i233.short