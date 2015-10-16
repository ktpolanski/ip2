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
	
#docker call, on kpolanski@tesla-wsbc.warwick.ac.uk
#-v for volume setup, the "real" directory has the proper path, then you get the colon
#and then you get the within-docker address (here /agave)
#-it so that it spits out the stdouts so we can see them
#and then just normal inputs. the working directory will be /agave, where the files are at
#--rm to kill the container once the thing completes to avoid littering the container space
docker run -it --rm -v /home/kpolanski/docker_demo/wigwams-image-building/testdata:/agave wigwams --Expression model_expr.csv --DEGs model_deg.csv --Export_Annotation annot_agi.tsv --Export_Hyperlink 'http://www.arabidopsis.org/servlets/TairObject?type=locus&name={gene}' --NoStandardising --PoolNumber 4 --SizeThresholds 10 10 8 5 5 --Job model --Legacy