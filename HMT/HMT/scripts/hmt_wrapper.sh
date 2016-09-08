#!/bin/bash
set -e

#inputs:
# $1 - logos folder for Paul's thing
# $2 onwards - python hypergeometric test arguments
# this includes $2 - a flag --Webapp if a webapp is to be made (if it's not to be made, then $2 will be something else)

#run the thing
python3 /scripts/hmt.py ${@:2}

#are we making a webapp
if [ $2 == '--Webapp' ]
then
	#and now time for Paul's logos thing
	cp -r $1 logos2 #no need to blank this later as it's moved within the script
	mkdir html
	perl /scripts/cyverse_hypergeometric_motif_test.pl . html

	#placeholder file cleanup to make iPlant output cleaner
	rm Webapp_*
fi