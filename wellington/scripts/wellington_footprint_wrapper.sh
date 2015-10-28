#!/bin/bash

#run wellington.
#note that we need to make the directory ourselves because yes
#we also need to run the indexing ourselves because yes
mkdir analysis
samtools index $2
python3 ../scripts/wellington_footprints.py $@

#move wellington run contents into main folder
mv analysis/* .
rmdir analysis

#get the average profile thing
mkdir output_visualisation
python3 ../scripts/dnase_average_profile.py WellingtonFootprints.FDR.bed $2 output_visualisation/average_footprint.png

#get the heatmap
python3 ../scripts/dnase_to_javatreeview.py WellingtonFootprints.FDR.bed $2 output_visualisation/javatreeview_heatmap_ready.csv

#get the wiggle tracks
python3 ../scripts/dnase_wig_tracks.py $1 $2 output_visualisation/fw_cuts.wig output_visualisation/rv_cuts.wig