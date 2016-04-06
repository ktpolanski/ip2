# Gaussian Process Two-Sample (GP2S) Test of Differential Expression

## The Purpose of the Algorithm

If available data takes on the form of a time course, with two separate conditions (such as control and treatment) present throughout the duration, then it would be preferable to be able to account for all of that information when evaluating differentially expressed genes. Additionally, such data could be mined for information on the timing of the divergence of the expression of the two conditions. GP2S is capable of performing both of those analyses.

## Basic Input/Output

GP2S accepts expression data as a CSV, with the first column being gene names, the first row delineating the two conditions (with the appropriate condition name appearing for each column of expression data), and the second row containing information on the time point. The expression data is assumed to be normally distributed, if you've got count data then log-transform your data before feeding it into the program.

The output comes in the form of a tab-delimited file with the gene names as the first column and the log Bayes factor as the second column. The higher the log Bayes factor, the more evidence there is for the gene being differentially expressed [(Kass and Raftery, 1995)][bayesfactor], and setting a cutoff on it (usually somewhere between 3 and 10) can produce an explicit list of differentially expressed genes. If the time-local mode is in use, the score file will also contain the identified time of first differential expression, and a supplemental file will be generated containing the probability of the gene being differentially expressed at each time evaluated by the algorithm.

## How Does It Work?

GP2S is based around the idea of fitting Gaussian processes to the time course data, producing two separate models to evaluate whether there is a difference in expression between the two conditions. One model features a single Gaussian process fit to all of the data at once, irrespective of condition, whilst the other fits each condition separately. The quality of the individual models is measured through the marginal likelihood, and subtracting the single-Gaussian-process log likelihood from the one-Gaussian-process-per-condition log likelihood yields the log Bayes factor, which is reported as the score to judge differential expression on.

The time-local extension sees the algorithm's idea applied in a different manner to try to elucidate the time at which differential expression begins to occur. The entirety of the experiment duration is densely divided up into equally spaced intervals, to get higher resolution than what the data itself brings to the table, and a Gibbs sampler is ran, optimising the allocation of each of the denser time scale points to either differentially expressed or not differentially expressed. At the end of the run, a probability is assigned to each of the points being differentially expressed, and the first point at which the probability crosses a provided threshold is reported as the time of first differential expression.

For further details, consult [Stegle et al., 2010][stegle2010].

## Test Run

If you want to get a feel for how GP2S operates without using your dataset, a small-scale demonstration dataset is provided at `ktpolanski/gp2s_testdata/input.csv` under Community Data. You can leave all the parameters as defaults, and check the time-local checkbox as appropriate if you want to test out the timing of first differential expression variant of the analysis.

## General Input

### Expression CSV

**Mandatory input.** The file containing your expression data that's to be mined for differential expression information. Gene names in the first column, first row for condition information (for example, control/treated), second row for time point information. Each condition needs to have the same time points, and each time point needs to have the same number of replicates across both conditions and all other time points. Refer to the demo file `ktpolanski/gp2s_testdata/input.csv` under Community Data for formatting if needed.

### Gaussian Process Hyperparameters

**Default:** 0.5,1,0.4

The starting values of the hyperparameters of the Gaussian processes used to fit the data, provided as a comma-delimited list of length scale, process variance and noise variance. These values are subsequently adjusted within the script to optimise the fit.

### Hyperparameter Optimisation Iterations

**Default:** 20

The hyperparameter optimisation is done through Expectation Propagation, with this parameter controlling the number of iterations of optimisation performed.

## Time-Local Input

**The following parameters will only be in use if GP2S is ran in time-local mode**

### Use Time-Local Mode?

If checked, the algorithm will perform an additional layer of analysis, aiming to elucidate the timing of differential expression of each gene. The standard test is also performed, with the log Bayes factor being reported as in the normal mode. However, the time-local mode is more computationally intensive than the standard test, which leads to this being an option and not a default.

### Prior Belief in Differential Expression

**Default:** 0.3

The prior belief that a gene is differentially expressed at any one time point, used to weight the individual likelihoods in the likelihood ratio to determine the reported probability of differential expression or lack thereof. At 0.5, the belief in differential expression is the same as the belief in the lack thereof and the weights cancel out. Higher values skew the test towards reporting differential expression.

### Time Point Resolution

**Default:** 20

The time-local version of GP2S takes the minimum and maximum time points reported in the data input and form an equally spaced, higher resolution time scale to evaluate differential expression on. This argument controls how many such evenly spaced time points should be generated.

### Initial Non-Differential Expression

**Default:** 2

Typically, time course experiments are designed in a manner that show the full temporal evolution of a given process. As such, it is likely that the start of the time course will have the two conditions behaving the same way. This parameter controls how many initial time points on the high-resolution, linearly spaced scale mentioned above are set to non-differentially-expressed when initialising the sampler designed to determine the probability of differential expression at each time point.

### Gibbs Sampler Iterations

**Default:** 30

The differential expression status of each time point on the higher resolution time scale is assessed with the help of a Gibbs sampler, with this parameter controlling the number of iterations of said sampler.

### Z Curve Hyperparameters

**Default:** 0.7,2,1E-5

Once the Gibbs sampler is done, the final probabilities of differential expression are smoothed with a second Gaussian process fit. As in the case of the Gaussian processes used to fit the data, the hyperparameters for this Gaussian process are accepted on input as a comma-delimited list of length scale, process variance and noise variance. Unlike the other Gaussian process hyperparameters, these are not tuned within the script.

### Differential Expression Probability Threshold

**Default:** 0.5

Upon smoothing the curve, a time of first differential expression can be elucidated by identifying the first point at which the probability of differential expression is at least this much. Increasing this parameter will lead to a more stringency on differential expression detection, whilst decreasing it will have the opposite effect.

## Output in Detail

### `scores.txt`

The primary output of the script, featuring the gene names as the first column and the log Bayes factor [(Kass and Raftery, 1995)][bayesfactor] scores as the second column. These can be used as a measure of differential expression, the higher the more confident GP2S is in the differential expression of the gene. Applying a cutoff (recommended somewhere between 3 and 10) will result in a simple, binary list of differential expression status. In the case of using the time-local mode, a third column will be added with the estimated time of first differential expression provided where applicable, and NaN where none was identified.

### `Z.txt`

Output specific to the time local mode, featuring the complete sequence of probabilities of differential expression at every point of the high-resolution, linearly spaced time scale. The time scale itself is not provided due to stability issues, but can be obtained easily given the fact it's a number of linearly spaced time points starting at the experimental start time and ending at the experimental end time time.

### `plots/`

A directory featuring a plot per gene, with the gene name being the plot file name. The blue fit is the joint fit across all the data irrespective of condition, while green and red are the individual condition fits. The solid line is the Gaussian process mean, while the shaded area around it has a width of two standard deviations. The log Bayes factor is also present on the plot. In the case of the time-local mode, the fits will often exhibit drastic shifts in the width of the shaded areas. This stems from the fact that during sampling, portions of the data are allocated to the single-Gaussian-process fit and individual-Gaussian-processes fits separately, and the regions allocated to one are not seen by the other, making it difficult for the fit to approximate what is happening in sections without allocated data. The probability of differential expression, as exported in `Z.txt`, is plotted in an additional subplot at the top of the time-local plots.

[bayesfactor]: http://amstat.tandfonline.com/doi/abs/10.1080/01621459.1995.10476572
[stegle2010]: http://online.liebertpub.com/doi/abs/10.1089/cmb.2009.0175