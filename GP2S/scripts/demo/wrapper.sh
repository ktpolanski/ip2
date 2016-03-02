#!/bin/bash

#IMPORTANT!
#when making the iPlant app, have the input CSV be the first thing in the argument order
#and have it be positional, i.e. no flag

#the GP2S script needs to run in its magical home directory or else everything fails because paths
SCRATCH=$PWD
cp $1 /scripts/demo/$1
cd /scripts/demo
mkdir plots

#anti-crash loop. gotta love randomly unstable code
until python run_two_sample.py $@; do
   python CSV_fix.py $1
   sleep 1
done

#move the output to our original working directory so iPlant can lap it up
cd $SCRATCH
if [ -f /scripts/demo/Z.txt ]; then
   mv /scripts/demo/Z.txt Z.txt
fi
mv /scripts/demo/scores.txt scores.txt
mv /scripts/demo/plots plots