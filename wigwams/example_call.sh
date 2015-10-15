#this runs the legacy version of Wigwams, as in the Polanski et al 2014 paper
#for slight improvements to the mining/merging, remove the --Legacy at the end
#for less redundant output, at the cost of some information, add --NonRedundantOutput
python wigwams_wrapper.py \
	--Expression model_expr.csv \
	--DEGs model_deg.csv \
	--Export_Annotation annot_agi.tsv \
	--Export_Hyperlink 'http://www.arabidopsis.org/servlets/TairObject?type=locus&name={gene}' \
	--NoStandardising \
	--PoolNumber 4 \
	--SizeThresholds 10 10 8 5 5 \
	--Job model \
	--Legacy