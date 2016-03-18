#!/bin/bash

#IMPORTANT!
#the first argument you need to provide is the expression CSV

cp $1 disposabledata.csv
mkdir plots

#anti-crash loop. gotta love randomly unstable code
until python run_two_sample.py disposabledata.csv ${2:@}; do
   python CSV_fix.py disposabledata.csv
   sleep 1
done