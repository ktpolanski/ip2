# hCSI

## The Purpose of the Algorithm

hCSI (hierarchical Causal Structure Inference) is an algorithm designed for inferring regulatory network models from time course expression datasets showing the same organism responding to different treatments/perturbations. The underlying assumption is that the response across the conditions may show some degree of similarity, with elements of a common, shared regulatory network (referred to as a hypernetwork) in use during the individual responses. The distinctness of the stimuli is tuned during the hCSI run, with data for conditions deemed dissimilar having less impact on the hypernetwork and, by extension, other conditions. The algorithm produces network models for each stimulus, as well as an estimate of the hypernetwork.

## Basic Input/Output

hCSI accepts normally distributed data on input, so if you want to use RNA-Seq data on input consider log-transforming it beforehand. The output of hCSI is a collection of CSV files with the edge probabilities in the network models for each perturbation, as well as the hypernetwork. These can be converted to GML files (which are compatible with a number of programs, including Cytoscape) with the hCSI_MarginalThreshold app.

## How Does It Work?

The basis for hCSI is CSI (Causal Structure Inference), a network inference algorithm designed to propose a regulatory network model for a single time course dataset. CSI operates under the premise that when a particular parent gene regulates a target (child) gene, then the child's expression profile can be explained by the parent's expression profile at the preceding time point (giving the parent gene time to get its protein to the child gene's promoter and affect its expression). This is formally done by fitting a Gaussian process to a visualisation of the data where one axis is made up of the child's expression and the other axes are time shifted parents. CSI infers parents on a per-child basis, with all eligible parent combinations (up to a specified maximum number of parents) evaluated, and probabilities of each parental set being "the correct one" are computed based on how good that individual fit is. The fits are optimised through Expectation Maximisation, with the values of hyperparameters governing the behaviour of the Gaussian processes tuned to maximise the overall fit quality. CSI was applied to synthetic 10 and 100 gene networks in a comparison with several other algorithms, and performed very well [(Penfold et al., 2011)][penfold2011].

hCSI takes CSI as a base and creates an individual CSI fit for each perturbation's worth of data fed into it. CSI would also be capable of analysing the data, but would produce a single network jointly inferred for all the data at once, whilst hCSI allows for a degree of deviation between the individual stimulus networks. This is achieved through the introduction of the hypernetwork, which can be interpreted as the underlying "master" network that the individual responses are based on, with an additional temperature parameter controlling the level of the hypernetwork's influence on a particular stimulus response. In contrast to CSI and its use of Expectation Maximisation, hCSI makes use of sampling, with a single step involving Gibbs sampling of individual dataset networks and the overarching hypernetwork, and Metropolis-Hastings sampling of the individual dataset hyperparameters and temperatures. Similarly to CSI, hCSI performs these procedures on a per-gene basis, with every gene taking on the role of the child whose expression is to be explained. The underlying mathematics can be found in [Penfold et al., 2012][penfold2012].

## Test Run

If you want to take hCSI out for a spin without using your own data, this can be done with the aid of one of the 10-gene synthetic networks originally used in the CSI and hCSI publications. The dataset to be used on input can be found at `ktpolanski/hcsi_testdata/dream4_5.csv` under Community Data. Leave all the parameter values as defaults. hCSI is computationally intensive and it may take upwards of an hour for this analysis to be complete.

## Input in Detail

**NOTE: hCSI is extremely computationally intensive.** It is greatly recommended to perform a high degree of preliminary analysis and select the most relevant subset of candidate genes to perform hCSI on. Going above 30 genes is really not recommended.

### Gene Expression CSV

**Obligatory input.** Comma-delimited file, with expression data ordered to have genes as rows and time points as columns. In terms of headers, the first column should contain gene IDs, the first row should contain condition names (repeated for each time point part of the condition), and the second row should contain the corresponding time of the time point in that condition. For reference on formatting, consult `ktpolanski/hcsi_testdata/dream4_5.csv` under Community Data.

### Parental Set Depth

**Default:** 2

When evaluating parental set combinations, a limitation is put on up to how many parents to sample from the parent pool to create the combinations. As this depth is increased, the number of parental sets to evaluate drastically goes up, making the problem less computationally tractable. Increasing this value above 2 is not recommended, especially if the dataset is larger than the test data provided. At the same time, lowering this value to 1 will result in uninformative analysis results, missing out on a lot of combinatorial regulatory action.

### Gaussian Process Prior

**Default:** 10,0.1

Part of CSI/hCSI is placing an assumption on how we expect the hyperparameter values to be distributed, and in both the publications the hyperparameters were assumed to be Gamma distributed with a shape parameter of 10 and a scale parameter of 0.1. If you wish to alter the shape/scale parameters of the gamma distribution, provide them as shape and scale with just a comma between them.

### Temperature Prior

**Default:** 1,1

The temperature parameters, introduced in hCSI, are also expected to be Gamma distributed, with the publication using a shape value of 1 and a scale value of 1. If you wish to alter it, the formatting is the same as for the Gaussian process priors.

### RNG Seed Offset

**Default:** 0

In order to avoid using the same exact random numbers for sampling when analysing different potential child genes, the RNG is initialised based on the row number of the child gene in the CSV expression file. This makes the output be deterministic for any individual data file provided on input. In case a change in the used RNG chains is desired (for example, in case the sampler gets stuck in a suboptimal local region), changing this parameter from 0 will lead to shifting the seeds by the provided value.

### Sample Count

**Default:** 25,000

The number of complete sampling cycles (Gibbs sampling of the individual networks and hypernetwork, Metropolis-Hastings sampling of hyperparameters and temperatures) to perform. Increasing this value will make the algorithm take longer to run, whilst decreasing it introduces the possibility that the algorithm won't manage to converge around the optimal solution.

### Burn-In Count

**Default:** 2,500

The algorithm is fully randomly initialised, so it will take it a while to arrive in the ballpark of the proper model. These many first samples will be discarded, writing them off as the time the sampler to land more less where it should be. Increasing this value will increase confidence in the convergence of the final result, but may require compensation from the sample count to ensure enough information is making it to the final model. Decreasing this value is not recommended.

### Data Normalisation

CSI fits are calibrated to expression data with zero mean and unit variance. By default, CSI standardises the expression data it receives on a per-gene, per-condition basis to match this requirement. If you believe you have strong reason to not standardise the data automatically within the code, you can choose to merely center it, or perform no normalisation whatsoever and leave it as is in the input.

### Gibbs Chain Storage

Samplers have their limits, and sometimes it's possible for them to get stuck in local optima wich may not actually be the best solution. Checking this box will produce a Python Pickle of the value chains of the individual networks and hypernetwork for possible diagnostic aid. The Pickled variable will be a list of lists, with a single item of the master list corresponding to an individual gene. These may be out of order relative to the input, as hCSI is parallelised for efficiency. An individual gene's list entry will in turn be a list of N+1 lists, where N is the number of conditions in the input file. The N+1 lists will be the conditions, sorted alphabetically, and the hypernetwork. Those lists will in turn feature all the parental set values with burn-in discarded. At this point it is possible to distinguish which gene the list of values is for, as it will always be the second element in the individual element tuple (with the first element being a list of parents). Storing the chains is not essential for every run, and will take up roughly 2-3 MB of space per gene in the dataset.

## Output in Detail

### `hcsi-XXXXX.csv`

One of these files will be generated per condition in the input dataset. This is the regulatory network model, with parents as columns and children as rows, and each edge given a probability score (the higher the better). Introduce a threshold (for example, in the hCSI_MarginalThreshold app) to turn this into a standard, binary network.

### `hcsi-hypernetwork.csv`

Exactly the same formatting as the individual condition files, features information on the identified hypernetwork. This might be of some use, so including it in the output - in the case of the 10-gene synthetic network provided as test data, the hypernetwork does a better job of reproducing the true network than any of the five individual models.



# hCSI_MarginalThreshold

## The Purpose of the Algorithm

hCSI_MarginalThreshold is a simple post-processing script which takes the marginal CSVs created by the hCSI app and turns them into GML files while applying a desired edge stringency threshold. The resulting GML files are compatible with a number of visualisation and analysis programs, including Cytoscape.

## Test Run

If you want to try out the app, but don't have access to any hCSI results of your own yet, feel free to use all the `ktpolanski/hcsi_testdata/hcsi_output/hcsi_XXX.csv` files under Community Data. These are the hCSI results for the demonstration `dream4_5.csv` dataset, with one file for each of the perturbations (A through E) and the hypernetwork. The computational load of this script is negligible in comparison to hCSI proper, so feel free to experiment with different values of the threshold and assess the resulting networks in any desired downstream software of choice.

## Input

**All input is obligatory, the app doesn't come with any defaults.**

### Marginal CSV Files

The CSV files that the hCSI app produces. Rows for children, columns for parents, first row and column as headers, each element representing the probability (from 0 to 1) that, according to hCSI, a particular parent regulates a particular child under the corresponding perturbation. It is fine to provide multiple files on input, such as a whole hCSI run's worth.

### Probability Threshold

The desired stringency for the GML network. Network edges below the provided probability threshold are removed from the export. The higher the threshold the fewer edges will be present in the GML network, but the higher hCSI's confidence in each one.

## Output

### `hcsi-XXXXX.csv.gml`

For each CSV file on input, a corresponding GML file will be created. These carry a rather minimal set of information - that the graph is directed, a list of nodes (using the gene names from the CSV file used on input for the underlying hCSI run) and a list of edges between them. The edges come with the marginal probability stored under the name of weight, and an extra parameter under the name of interaction is set to "pd" (protein-DNA) in the interest of Cytoscape friendliness.

[penfold2011]: http://rsfs.royalsocietypublishing.org/content/1/6/857.short
[penfold2012]: http://bioinformatics.oxfordjournals.org/content/28/12/i233.short