# Hypergeometric Motif Test

## The Purpose of the Algorithm

When a group of genes of interest is obtained, be it through computational (such as clustering of expression profiles) or experimental (such as ChIP-Seq) means, it is often of interest to assess its biological functionality. This can be done by taking known biological information, such as GO terms, signalling pathways, or transcription factor binding site presence in promoters, and assigning them to genes across the genome and subsequently evaluating the overrepresentation thereof in the group of interest. As an example, [BiNGO][bingo] is capable of evaluating GO term overrepresentation using the hypergeometric motif test. This set of two apps allows you to do the same, using known transcription binding sites in place of GO terms. The first app, `HMT_Index`, creates an index of motif hits using [FIMO][fimo] and statistical postprocessing, and is quite computationally intensive. However, it only needs to be ran once, as the resulting index file can be used on input for any number of `HMT` analyses for that particular organism and motif set, with the second script in the pipeline being far quicker.

# `HMT_Index`

## The Purpose of the Algorithm

When [BiNGO][bingo] performs its GO term overrepresentation analysis, it makes use of a comprehensive GO term annotation which links individual GO terms to particular genes. Armed with this information, it can quite easily compute the sizes of the two sets that are being scanned for an overlap (the number of genes with a particular GO term versus the size of the gene group of interest), identify said overlap, and perform a hypergeometric test to assess said overlap's significance. In the case of transcription factor binding sites, no such unified annotations for the presence or absence thereof in promoters of genes exist, in no small part thanks to the lack of a definite motif collection and variable desired promoter sizes.

The `HMT_Index` app creates the equivalent of a GO term annotation, using transcription factor binding motifs of interest provided on input. The program scans the promoters of genes for motif hits using [FIMO][fimo], and then performs a statistical follow-up analysis to assess whether a binding site is deemed to be present or absent from a particular promoter. The follow-up makes use of [the Hommel method][hommel] to combine up to five FIMO-identified hits' p-values, basing on the assumption that multiple weaker hits of a particular motif in a promoter can be as relevant as a single strong one, into a single p-value. Said p-value is used in a binomial test to determine whether the motif is deemed present or absent from the promoter of choice.

## Test Run

If you want to try out the indexing part of the hypergeometric motif test, the repeat-masked Arabidopsis genome FASTA file and TAIR10 GG3 annotation are provided in `ktpolanski/hmt_index/testdata` along with a single MYB52 binding motif [(Franco-Zorrilla et al., 2014)][fz] in MEME formatting.

## Input

### Genome

The app requires a genome sequence on input so that promoter sequences can be easily automatically extracted. Please provide the genome as a single FASTA file.

### GFF3 Annotation

The GFF3 annotation needs to be compatible with the provided genome, and is used to locate the genes for promoter extraction. In the second app, the GFF3 annotation is used to create a complete universe to use as the background in hypergeometric testing.

### GFF3 Gene ID Attribute

A GFF3 file can carry a lot of information about an organism's genes, whilst the program is after the very basics - a distinct and discernible gene ID style that is also used in the input of gene groups for the second app, along with corresponding positioning. The GFF3 file is filtered to the lines that contain gene information, but the script subsequently needs information on which of the information fields to use as the identifier. For example, the Arabidopsis test data provided (`ktpolanski/hmt_index_testdata/annot.gff3` under Community Data) has the `gene_id=` field correspond to AGI identifiers, which are the widely accepted locus code nomenclature for Arabidopsis.

### Promoter Length

How many base pairs upstream of the transcription start site to take as the promoter sequence for analysis. Some common values include 200, 500 and 1000 base pairs.

### Gene Overlap Removal

If genes are close together, it is possible that the designated promoter region for one gene is going to overlap with the transcribed area of another gene. This does not necessarily disqualify regulatory motifs being present in the overlapping region. Checking this box will trim the parts of the promoter regions that are also part of the transcribed areas of other genes.

### Motif Input

The transcription factor binding site motifs you wish to identify in the genes' promoters. These need to be either MEME input or Uniprobe formatted. Demonstration databases can be downloaded from [the MEME Suite page][meme]. The page also offers a number of utility scripts to convert motif sequences from other formats into MEME input.

### Uniprobe Formatting

Uniprobe support is present in the app as Uniprobe is, arguably, the simplest and most intuitive format of motif formatting. This makes it fit for manual motif construction or easy reformatting if so desired. Checking this checkbox signifies that the input is Uniprobe formatted, and the script will automatically convert it into MEME input.

### Significance Threshold

The desired significance cutoff point. P-values below this threshold are deemed significant, first during FIMO motif mining, then during the statistical post-processing. The second stage features Bonferroni false discovery rate correction to ensure high stringency in motif detection. 0.05 is widely accepted and used.

### Top Motif Count

As previously mentioned, an assumption of this script is that many weaker motif hits may be just as effective as a single stronger hit. This stems from transcription factors' propensity to become detached from their binding sites and slide along the DNA, with the possibility of getting bound again at another binding site. The top non-overlapping motif hits will be sorted on p-value and up to this many p-values will be combined into a single p-value for use in the binomial follow-up test. Arabidopsis testing has fared well with taking the top 5 significant hits per motif per gene.

## Output

### `fimo_found.txt`

The app produces a single file, which features two columns and is tab delimited. The first column is a motif name matching one of the motifs provided in the input file, whilst the second column is a gene ID. For any given motif-gene pair in the output, the statistical post-processing of the FIMO output found some combination of the top hits of the motif in the gene's promoter to be significant. This file is to be used on input for the `HMT` app.

### `logos/`

A folder including PNG images for the PWMs of each of the motifs used in the mining. Taken on input by `HMT` to produce an interactive visualisation webapp.

### `universe.txt`

A file containing all the gene IDs within the analysed organism, taken on input by `HMT` and used as part of the statistical evaluation of motif overrepresentation. Contains all the gene IDs within the tested organism, with one gene ID per line.



# `HMT`

## The Purpose of the Algorithm

The second app in the pipeline, `HMT` begins where `HMT_Index` ends and uses the produced "motif annotation" to perform overrepresentation analyses of gene groups of interest. The same `fimo_found.txt` file produced by `HMT_Index` can be used on input for as many analyses for the same organism and motif set as desired.

The algorithm itself is a very straightforward hypergeometric test - given a particular universe and two sets (genes in the gene group versus genes with a particular motif), an overlap is identified and its significance is assessed. The returned output files are very comprehensive, featuring a number of different formats and FDR corrections.

## Test Run

If you want to get a feel for the output of the hypergeometric motif test tool, you'll find all the required files at `ktpolanski/hmt_testdata` under Community Data. The `fimo_found_200bp.txt` file is a motif hit file generated using the `HMT_index` app for a collection of multiple protein binding microarray-derived motifs for Arabidopsis transcription factors ([Franco-Zorrilla et al., 2014][], [Weirauch et al., 2014][weirauch]), while `input.txt` is a collection of 78 gene groups.

## Input

### Gene Group Input File

A tab-delimited file with two columns, with the first column being the number of the gene group and the second column being a single gene ID in that gene group. For example, if you wanted to analyse two gene groups, with the first one containing 10 genes and the second one containing 15 genes, the input file would be 25 lines long, with 10 lines of gene group 1 IDs and 15 lines of gene group 2 IDs. Consult an example Arabidopsis file at `ktpolanski/hmt_testdata/input.txt` under Community Data.

### Motif Hit File

`fimo_found.txt`, as produced by `HMT_Index`. You can also supply your own motif hit file if so desired - the first column is to feature a motif name with the second column featuring the gene ID where the motif was sighted, with a tab delimiting them. An example Arabidopsis file can be found at `ktpolanski/testdata/fimo_found.txt` under Community Data.

### PWM Logos Folder

The `logos` folder, as produced by `HMT_Index`. If you provide your own motifs, have this folder include a PNG PWM image for each motif, named with the motif named. An example Arabidopsis folder can be found at `ktpolanski/testdata/logos/` under Community Data.

### Gene Universe

The `universe.txt` file produced by `HMT_Index`. If you're supplying your own gene hit files, then create the universe as a collection of all the gene IDs in the organism, with one line per gene ID. An example Arabidopsis file can be found at `ktpolanski/testdata/universe.txt` under Community Data.

### Significance Threshold

The p-value cutoff point for a result to be deemed significant. A complete listing of p-values, both significant and not, is returned in case it's of use.

### Make Webapp

If the box is checked, an interactive web browser of the raw p-value results will be produced. In the case of large searches, this will greatly slow down the generation of the results due to extensive parsing, with the standard outputs providing a focused list of relevant overrepresentations.

### Motif Annotation

If desired, the output can be enhanced with extra information on the motifs used in the search, such as which transcription factors bind said motif. Optional input. If provided, it's to be a CSV file with no header (row/column names), with the first column being the motif ID that is to match the ones found in `fimo_found.txt`. The annotation will be ignored if any duplicate entries are identified or any of the `fimo_found.txt` motif IDs are not found, but it will not crash the app.

## Output

The resulting output files are created after three separate FDR corrections, allowing the user flexibility in terms of which one to use in later analyses:
* **Benjamini-Hochberg** - a Benjamini-Hochberg correction is performed on a per-gene-group basis, the least stringent of the featured corrections
* **Bonferroni** - a Bonferroni correction is performed on a per-gene-group basis
* **GlobalBonferroni** - a single Bonferroni correction is applied across all gene groups, the most stringent of the featured corrections

### `Full_XXX_P-Values.txt`

A complete P-value matrix, with rows as motifs and columns as gene groups. Features both the significant and insignificant test results. Provided for all three FDR corrections, along with a matrix of raw p-values as well.

### `Significant_XXX_P-Values.txt`

The previous file group, but filtered to feature the rows and columns where at least one item is significant. This compressed version is good for heatmap creation (for example, in Excel with the aid of conditional cell formatting). Provided for all three FDR corrections.

### `Overrepresentation_XXX.txt`

A comprehensive listing of all the significant overrepresentation instances, featuring information on the motif ID, gene group ID, corrected p-value and a complete list of all the genes from the gene group that had the motif identified in their promoters. Provided for all three FDR corrections.

### `html/`

A folder containing an interactive webapp allowing the browsing of raw p-values. The cells can be shaded based on p-value thresholds, and individual overrepresentation p-values are clickable to reveal a list of genes in the gene group that had a given motif detected in their promoters. Only produced if the `Make Webapp` box is checked.

### `FullOutput.tar`

The complete output of the analysis, archived into a single file for ease of downloading to your computer.

[bingo]: https://bioinformatics.oxfordjournals.org/content/21/16/3448.full
[fimo]: https://bioinformatics.oxfordjournals.org/content/27/7/1017.full
[hommel]: http://arxiv.org/pdf/1212.4966.pdf
[meme]: http://meme-suite.org/db/motifs
[fz]: http://www.pnas.org/content/111/6/2367.short
[weirauch]: http://www.cell.com/cell/abstract/S0092-8674(14)01036-8?_returnURL=http%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867414010368%3Fshowall%3Dtrue&cc=y=