# Wigwams Identifies Genes Working Across Multiple Situations

## The Purpose of the Algorithm

Wigwams could be described as a biclustering approach for time series data, as it identifies genes exhibiting co-expression across subsets of available time course data sets. This is a sensible course of action when looking for co-regulatory crosstalk in the responses to multiple stimuli, as just because you feed an algorithm data from a high number of experiments does not oblige the biological crosstalk to span across all of them. 

In contrast to standard biclustering methods, Wigwams evaluates the co-expression trends it identifies for statistical significance, discriminating between co-expression indicative of co-regulation and co-expression occurring by chance, merely stemming from the relative abundance of profiles in each individual condition. Once again, the biological motivation is sound - if given 5000 genes exhibiting a particular profile in one condition and 5000 genes exhibiting another particular profile in a second condition, an overlap of 50 between the two sets of 5000 does not seem indicative of crosstalk and is probably occurring by chance. By contrast, if the individual condition profiles occurred 100 times instead of 5000, the 50 gene intersection would be considerably more likely to represent a shared regulatory mechanism. Existing biclustering methodology does not attempt to discriminate between the two scenarios, and will likely identify and report the 50 gene group without discriminating on the context of the co-expression.

Additionally, Wigwams is capable of accounting of genes' differential expression status across the conditions, making sure that the crosstalk it detects is performed by genes playing an active part in the responses.

## The Components of a Wigwams Run

When performing a Wigwams analysis, the data is first mined for groups of genes, henceforth termed modules, exhibiting co-expression indicative of co-regulation across subsets of the provided time course datasets. The resulting module list features a high degree of redundant information, so the following step is to perform post-processing to make the output be succinct and accessible.

During module mining, all of the genes differentially expressed in at least two of the supplied experiments take turns being the "seed gene" - a reference profile that the expression of all the other genes is compared to. For each individual condition where the "seed gene" is differentially expressed, a list of the genes exhibiting the strongest co-expression with the "seed gene" profile is compiled. The individual lists are then compared to each other in all possible combinations, with the overlap being evaluated for significance with a modification of the hypergeometric test (a tool used for general overrepresentation assessment, for details consult [Polanski et al. 2014][polanski2014]). Stringent FDR correction is performed by applying the Bonferroni correction to the resulting p-values.

As mentioned, the module mining evaluates all of the eligible genes as "seed genes", with all of the possible condition combinations mined for modules as well. This leads to a high degree of redundancy in the resulting module list, which needs to be dealt with to produce a clear output. Going back to the 50 gene module example, the complete search will find the module when using the first gene in the module as a "seed gene", then find it again when the second gene is used as the "seed gene", and so on in this fashion, resulting in many entries in the module list saying essentially the same thing. If the 50-gene module spanned three conditions instead of two, then the identified module list would also include the module detected across all corresponding pairs of conditions as well. At the time of module detection, the algorithm has no way of knowing that the three-condition module exists, and identifies everything that it finds evidence of co-regulation for. This is why we need merging, which glues together similar modules spanning the same conditions, and sweeping, which kicks out modules that are saying similar things to other modules spanning more conditions.

## Demonstration Run

You can quite easily reproduce the Wigwams run performed in [Polanski et al. (2014)][polanski2014]. The data can be found in `ktpolanski/wigwams_testdata` under Community Data. In the iPlant app, in the Inputs tab, use `model_expr.csv` for the `Gene Expression CSV` field, and `model_deg.csv` for the `Differential Expression CSV` field. In the Parameters tab, check the `No Standardising` and `Run Legacy Version` checkboxes, and enter `10;10;8;5;5` into the `Size Thresholding` field. Press the Launch Analysis button in the bottom right corner of the app window and your demonstration run has been launched.

## Inputs

### Gene Expression CSV

**The only obligatory file you have to provide.** After all, performing an analysis without any data is not really possible. This is where your data goes, formatted into CSV (the individual fields being separated by commas). The first column is to contain the gene IDs, which should be the same as the gene IDs provided in the differential expression file if you choose to use one. The first two rows should describe the particular time point that the data is for, with the name of the condition in the first row and the actual time point in the second row. If in need of an example, consult `ktpolanski/wigwams_testdata/model_expr.csv` under Community Data for formatting. If in possession of multiple replicates, it would be highly preferable to average them out to a single mean value per gene per condition.

### Differential Expression CSV

**Default:** When no file is provided, everything is differentially expressed everywhere, i.e. no differential expression information taken into account during analysis.

A CSV capturing the differential expression status of your genes across the conditions. The first column is to contain the gene IDs, which should match the ones from the expression CSV you just provided. Order does not matter though. The first row is to feature the names of the conditions, which should match the ones that you provided data for in the previous step as well. The actual data fields are to be 0 if the gene is not differentially expressed in that condition, or 1 if the gene is differentially expressed in that condition. If in need of an example, consult `ktpolanski/wigwams_testdata/model_deg.csv` under Community Data for formatting. You actually have to perform the differential expression analyses yourself, though. If in need of method guidance, GP2S ([Stegle et al., 2010][stegle2010]) is a solid algorithm, and it was used to identify differentially expressed genes in five of the six datasets comprising the demonstration files (the sixth dataset lacked a control time course to compare the treatment to, as the experiment was plant ageing).

### Gene Annotation

**Default:** No annotation. All exports are made using the IDs provided in the input CSV(s).

When looking at the exported membership of each of the modules, merely having the IDs that you provided in the CSV on input may not be the most immediately informative. As such, you can provide extra information on each of your genes to include in the export to help make sense of the modules faster. If you choose to provide an annotation, it needs to be tab-separated, with the first column matching the IDs that you provided in your CSV, the second column being a widely recognised form of gene ID (for example, your CSV could feature microarray probe names as IDs, and this column would feature actual gene identifiers), with columns three onwards featuring any additional information you may wish to provide. **Column two needs to contain a widely recognised form of gene ID** - if you already feature that as your CSV ID, just copy it over again as column two. For reference on formatting, including the duplication of informative CSV gene IDs as column two, consult `ktpolanski/wigwams_testdata/annot_agi.tsv` under Community Data. **Do not include any unmapped IDs matching your CSV - if the script fails to find the ID of a gene, it will generate information on it being unmapped by itself. Including information on unmapped IDs messes with the BiNGO/MEME-friendly output.**

### Hyperlink

**Default:** No hyperlink.

If you provide an annotation, Wigwams can automatically create hyperlinks to online resources for you to quickly and easily examine any gene that piques your interest. When providing a hyperlink, the place where the gene ID is to go should be marked with `{gene}`. As such, an example hyperlink, which works for the *A. thaliana* resource TAIR, is `http://www.arabidopsis.org/servlets/TairObject?type=locus&name={gene}`. If our gene is AT1G25550, then the hyperlink address that will show up in the module export will be `http://www.arabidopsis.org/servlets/TairObject?type=locus&name=AT1G25550`. The hyperlink column will be inserted as the third column in the annotation. Open the resulting `exported_modules.tsv` file in Excel for the hyperlinks to become active.

## Parameters

### Mining - Set Sizes

**Default:** 50;100;150;200;250

When performing module mining, Wigwams evaluates the top co-expressed genes with a single "seed gene" to detect co-regulatory phenomena. However, such phenomena can be occurring at different scales, with the potential to lose some of the smaller scale ones (a significant overlap of a low number of top co-expressed genes, but no additional overlap as the number of evaluated co-expressed genes becomes increased). As such, Wigwams evaluates a variety of sizes of top co-expressed gene sets to avoid such situations. **If performing Wigwams on a small scale dataset, make sure that your largest set size does not exceed the total number of genes in your dataset.** Scale them back accordingly, and try to keep them relatively small in comparison to your overall dataset size (try not to exceed 50% of your gene pool, even in very small datasets). Provide your desired set sizes as a semicolon-delimited list.

### Mining - Alpha

**Default:** 0.05

The significance threshold for judging whether the expression you're analysing is statistically significant or not. Within the program, this is corrected using the Bonferroni method to ensure maximum stringency. Increasing this value will result in more lenient initial module detection, whilst decreasing the value will increase the stringency of the procedure.

### Mining - Correlation Net

**Default:** 0.7

Sometimes you can bump into a small scale regulatory phenomenon, where only a small number of genes exhibit a given behaviour. If this is the case, then it is quite probable that some of the top co-expressed genes captured in the pre-determined set sizes will not actually be that regulated with the "seed gene" anymore. As such, the correlation net is there as a precaution to avoid creating modules featuring genes not correlated with the "seed gene" profile - even if the not-quite-correlated genes are co-regulated, it's probably a different co-regulatory event than the one being detected at this particular "seed gene" and best left for another iteration of the search to discover. If the top co-expressed genes with the "seed gene" stop being at least this co-expressed with the "seed gene" (as measured in Pearson's Correlation Coefficient) across each condition, the top co-expressed genes identified for the corresponding set size will not be mined for modules. Increasing this parameter increases the stringency of this condition, whilst decreasing it allows less tightly co-expressed initial modules to be identified. In an ideal world, a correlation net-style parameter would control the entirety of the module mining, by demanding that genes be at least that much co-expressed with the "seed gene" profile, but due to the implementation of the statistical testing the variable set sizes this would result in would render the procedure extremely computationally intensive.

### Merging - Overlap

**Default:** 0.3

The mining's complete search approach results in a high degree of redundancy among the module list, and merging is a procedure to deal with redundancy among modules spanning the same time course dataset subset. All modules for the same time course dataset subset span are compared pairwise, and if at least this proportion of the smaller module's genes are also present in the larger module the pair of modules is deemed redundant and merging commences. Increasing the parameter results in more distinct modules per condition combination in the final output, with possibility for a greater degree of left over redundancy. Decreasing the parameter results in the final modules capturing broader regulatory phenomena and potentially losing tightness of the expression profiles within each condition.

### Merging - Mean Correlation

**Default:** 0.9

Sometimes a regulatory phenomenon is of a very large scale, and modules exhibiting very similar expression profiles do not get identified as redundant due to the sheer volume of the genes involved in the process. If this occurs, the lack of redundancy is more likely to stem from subtle noise in the data than actual differing regulatory processes. As such, merging can also trigger and identify a pair of modules spanning the same conditions as redundant if their mean expression profiles across all the conditions are at least this correlated (measured through Pearson's Correlation Coefficient). Increasing this parameter can result in the presence of similar-looking, albeit not heavily gene redundant modules in the output that are likely part of the same regulatory process, but with their expression slightly differentiated by measurement noise. Decreasing this parameter can result in the merging of modules stemming from slightly differing regulatory phenomena.

### Merging - Correlation Filter

**Default:** 0.8

When merging a pair of modules, the larger module is taken as the reference and genes from the smaller module are transferred over if they are sufficiently correlated (measured through Pearson's Correlation Coefficient) with the mean of the larger module in each condition. This is in place to try to control the rate at which the expression of the merged module loses tightness. Increasing this value will result in very stringent gene selection and information loss (measured as unique genes present in the output), but lead to preserved tight regulatory phenomena. Decreasing this value will result in information preservation, but less tight regulatory phenomena.

### Sweeping - Overlap

**Default:** 0.5

Whilst merging deals with redundancy among modules spanning the same time course dataset subsets, it leaves redundancy among modules spanning different time course dataset subsets. As a practical example, if a particular phenomenon spans three conditions, it will get picked up as a three condition module, as well as every subset pair of conditions. Let's name the module spanning more conditions Module A, and the module spanning a subset of Module A's conditions Module B. Module B will be removed if the proportion of the genes it shares with Module A crosses this parameter value. Increasing this value will result in inter-condition redundancy in the output, but will preserve potentially broader local regulatory phenomena. Decreasing this value will result in potential information loss, but a less redundant output.

### Size Thresholding

**Default:** None. We don't know how many time course datasets you'll be using. For the analysis in [Polanski et al. (2014)][polanski2014], we used 10;10;8;5;5

At the end of the Wigwams analysis, you may wish to remove small modules to not clutter the output with very petite phenomena. This is an optional step, and leaving the field empty will not filter your module list in any way. If you do wish to perform this filtering step, provide N-1 desired module sizes, where N is the total number of time course datasets used on input, with the desired module sizes going from 2 to N condition span. Provide your desired size thresholds as a semicolon-delimited list.

## Algorithm Alterations and Enhancements

### No Standardising

By default, Wigwams will scale the data to a zero-mean, unit-variance normal distribution on a per-gene, per-condition basis, but if you desire to keep your data in the form that is provided in the CSV then check this box. Disclaimer: Wigwams uses Pearson's Correlation Coefficient for all co-expression evaluation, so if your data is not centered (mean of zero on a per-gene, per-condition basis) then your export plots can get messy and obfuscate expression trends.

### Run Legacy Version

Since the publication of the algorithm in [Polanski et al. (2014)][polanski2014], slight improvements have been made to the module mining and merging procedures. In the original publication, the most significant module was identified from each set size combination (for a given "seed gene" and condition combination), but more information is present in the output if the largest significant module is returned instead. The original design of merging did not foresee the creation of redundant modules that would have satisfied merging criteria through independent merging procedures, so merging is repeated until no new merging instances occur. In order to run the algorithm as originally outlined in the Bioinformatics article, check this box.

### Non-Redundant Output

If producing very succinct and focused output is desired, the redundancy removal procedures can be employed in an alternate manner. If given a regulatory phenomenon spanning three conditions, the corresponding two-condition modules can be made less redundant if the full three condition module is reconstructed first and is used to remove un-merged, initial, small modules across the two-condition subsets. If this box is checked, the algorithm will iteratively perform merging, optional filtering and sweeping using modules spanning fewer and fewer conditions to remove small initial modules spanning fewer conditions to make the output less redundant. **This introduces an additional parameter into the mix, as outlined below.**

### Non-Redundant Output - Sweeping Lower Bound

**Default:** 0.1. Inactive outside of Non-Redundant Output mode.

In Non-Redundant Output mode, the entirety of the redundancy removal procedure is repeated using different sweeping thresholds, minimising redundancy of the output (measured as the ratio of number of unique genes present in the final modules to the total module size). This requires the specification of a sweeping ratio range to evaluate, with the original `Sweeping - Overlap` parameter serving as the upper bound. The resolution is 0.05. Tuning this parameter is a trade-off between lack of redundancy at the cost of information loss (lower values) and more information at the cost of the presence of redundancy (higher values).

## Output

### `exported_modules-analysis.tsv`

This is the tab-separated text file containing the gene membership for the modules that were identified. You can open the file in Excel. The first column contains the Module ID, which is a unique number assigned to the module. The second column features the condition span, the time course experiments across which the module members are dependently co-expressed. The third column features the ID of the gene (as in your CSV input) that's part of the module, and you get one line per gene in the module. If you supplied an annotation, the corresponding annotation columns (where applicable) will constitute columns four onwards. If you supplied a hyperlink template, those will be present as column five.

### `plots-analysis/`

Subfolder with one plot per module. The plot name is the module ID, matching the first column of `exported_modules-analysis.tsv`. The plots are vector graphics (.eps format), making them publication-ready. The conditions highlighted in red are the ones where the module is deemed co-regulated, with the black ones not being identified as part of the module's condition span.

### `functional_analysis_inputs-analysis/`

The module members, handily formatted for follow-up analysis in BiNGO and MEME-LaB. BiNGO ([Maere et al. 2005][maere2005]) is a Cytoscape plugin that can simultaneously mine multiple groups of genes for GO term overexpression, but requires specific data formatting. The MEME-LaB-friendly input can be used for overrepresented motif mining in the web tool MEME-LaB ([Brown et al. 2013][brown2013]). If you provided an annotation, the export will be based on the fourth column of `exported_modules-analysis.tsv` with any unmapped gene skipped. In the case of no annotation, the export will be based on the third column instead.

### `analysis/`

This folder houses the intermediate module processing stages, if those are of interest to advanced users. `raw_modules.tsv` features the initial, heavily redundant mined module list, `merged_modules.tsv` features the module list after merging, `swept_modules.tsv` features the module list after sweeping, and `filtered_modules.tsv` features the module list after filtering (if you performed that step in your analysis). The files are tab separated, with the condition span of the module in the first column, the "seed gene" in the second column, the set size in the third column, the log10 p-value in the fourth column, and comma-separated CSV ID gene membership in the fifth column. Columns two through four only carry information in `raw_modules.tsv`, and are only kept in the other files for formatting consistency. In the case of Non-Redundant Output mode, the files are from the final pass of the iterative redundancy removal.

[brown2013]: http://bioinformatics.oxfordjournals.org/content/29/13/1696.long
[maere2005]: http://bioinformatics.oxfordjournals.org/content/21/16/3448.long
[polanski2014]: http://bioinformatics.oxfordjournals.org/content/30/7/962.long
[stegle2010]: http://online.liebertpub.com/doi/abs/10.1089/cmb.2009.0175