# Wigwams Identifies Genes Working Across Multiple Situations

## The Purpose of the Algorithm

Wigwams could be described as a biclustering approach for time series data, as it identifies genes exhibiting co-expression across subsets of available time course data sets. This is a sensible course of action when looking for co-regulatory crosstalk in the responses to multiple stimuli, as just because you feed an algorithm data from a high number of experiments does not oblige the biological crosstalk to span across all of them. 

In contrast to standard biclustering methods, Wigwams evaluates the co-expression trends it identifies for statistical significance, discriminating between co-expression indicative of co-regulation and co-expression occurring by chance, merely stemming from the relative abundance of profiles in each individual condition. Once again, the biological motivation is sound - if given 5000 genes exhibiting a particular profile in one condition and 5000 genes exhibiting another particular profile in a second condition, an overlap of 50 between the two sets of 5000 does not seem indicative of crosstalk and is probably occurring by chance. By contrast, if the individual condition profiles occurred 100 times instead of 5000, the 50 gene intersection would be considerably more likely to represent a shared regulatory mechanism. Existing biclustering methodology does not attempt to discriminate between the two scenarios, and will likely identify and report the 50 gene group without discriminating on the context of the co-expression.

Additionally, Wigwams is capable of accounting of genes' differential expression status across the conditions, making sure that the crosstalk it detects is performed by genes playing an active part in the responses.

## The Components of a Wigwams Run

When performing a Wigwams analysis, the data is first mined for groups of genes, henceforth termed modules, exhibiting co-expression indicative of co-regulation across subsets of the provided time course datasets. The resulting module list features a high degree of redundant information, so the following step is to perform post-processing to make the output be succinct and accessible.

During module mining, all of the genes differentially expressed in at least two of the supplied experiments take turns being the "seed gene" - a reference profile that the expression of all the other genes is compared to. For each individual condition where the "seed gene" is differentially expressed, a list of the genes exhibiting the strongest co-expression with the "seed gene" profile is compiled. The individual lists are then compared to each other in all possible combinations, with the overlap being evaluated for significance with a modification of the hypergeometric test (a tool used for general overrepresentation assessment, for details consult Polanski et al. 2014). Stringent FDR correction is performed by applying the Bonferroni correction to the resulting p-values.

As mentioned, the module mining evaluates all of the eligible genes as "seed genes", with all of the possible condition combinations mined for modules as well. This leads to a high degree of redundancy in the resulting module list, which needs to be dealt with to produce a clear output. Going back to the 50 gene module example, the complete search will find the module when using the first gene in the module as a "seed gene", then find it again when the second gene is used as the "seed gene", and so on in this fashion, resulting in many entries in the module list saying essentially the same thing. If the 50-gene module spanned three conditions instead of two, then the identified module list would also include the module detected across all corresponding pairs of conditions as well. At the time of module detection, the algorithm has no way of knowing that the three-condition module exists, and identifies everything that it finds evidence of co-regulation for. This is why we need merging, which glues together similar modules spanning the same conditions, and sweeping, which kicks out modules that are saying similar things to other modules spanning more conditions.

## Demonstration Run

You can quite easily reproduce the Wigwams run performed in Polanski et al. (2014). The data can be found in `ktpolanski/wigwams_testdata` under Community Data. In the iPlant app, in the Inputs tab, use `model_expr.csv` for the `Gene Expression CSV` field, and `model_deg.csv` for the `Differential Expression CSV` field. In the Parameters tab, check the `Run Legacy Version` checkbox, and enter `10 10 8 5 5` into the `Minimum Module Size Filtering` field. Press the Launch Analysis button in the bottom right corner of the app window and your demonstration run has been launched.

## Inputs

### Gene Expression CSV

**The only obligatory file you have to provide.** After all, performing an analysis without any data is not really possible. This is where your data goes, formatted into CSV (the individual fields being separated by commas). The first column is to contain the gene IDs, which should be the same as the gene IDs provided in the differential expression file if you choose to use one. The first two rows should describe the particular time point that the data is for, with the name of the condition in the first row and the actual time point in the second row. If in need of an example, consult `ktpolanski/wigwams_testdata/model_expr.csv` under Community Data for formatting. If in possession of multiple replicates, it would be highly preferable to average them out to a single mean value per gene per condition.

### Differential Expression CSV

**Default:** When no file is provided, everything is differentially expressed everywhere, i.e. no differential expression information taken into account during analysis.

A CSV capturing the differential expression status of your genes across the conditions. The first column is to contain the gene IDs, which should match the ones from the expression CSV you just provided. Order does not matter though. The first row is to feature the names of the conditions, which should match the ones that you provided data for in the previous step as well. The actual data fields are to be 0 if the gene is not differentially expressed in that condition, or 1 if the gene is differentially expressed in that condition. If in need of an example, consult `ktpolanski/wigwams_testdata/model_deg.csv` under Community Data for formatting. You actually have to perform the differential expression analyses yourself, though. If in need of method guidance, GP2S (Stegle et al., 2010) is a solid algorithm, and it was used to identify differentially expressed genes in five of the six datasets comprising the demonstration files (the sixth dataset lacked a control time course to compare the treatment to, as the experiment was plant ageing).

### Annotation

**Default:** No annotation. All exports are made using the IDs provided in the input CSV(s).

When looking at the exported membership of each of the modules, merely having the IDs that you provided in the CSV on input may not be the most immediately informative. As such, you can provide extra information on each of your genes to include in the export to help make sense of the modules faster. If you choose to provide an annotation, it needs to be tab-separated, with the first column matching the IDs that you provided in your CSV, the second column being a widely recognised form of gene ID (for example, your CSV could feature microarray probe names as IDs, and this column would feature actual gene identifiers), with columns three onwards featuring any additional information you may wish to provide. **Column two needs to contain a widely recognised form of gene ID** - if you already feature that as your CSV ID, just copy it over again as column two. For reference on formatting, including the duplication of informative CSV gene IDs as column two, consult `ktpolanski/wigwams_testdata/annot_agi.tsv`.

### Hyperlink

**Default:** No hyperlink.

If you provide an annotation, Wigwams can automatically create hyperlinks to online resources for you to quickly and easily examine any gene that piques your interest. When providing a hyperlink, the place where the gene ID is to go should be marked with `{gene}`. As such, an example hyperlink, which works for the *A. thaliana* resource TAIR, is `http://www.arabidopsis.org/servlets/TairObject?type=locus&name={gene}`. If our gene is AT1G25550, then the hyperlink address that will show up in the module export will be `http://www.arabidopsis.org/servlets/TairObject?type=locus&name=AT1G25550`. The hyperlink column will be inserted as the third column in the annotation.