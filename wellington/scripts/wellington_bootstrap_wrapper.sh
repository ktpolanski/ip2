#!/bin/bash

#run wellington.
#we need to run the indexing ourselves because yes
samtools index $1
samtools index $2
#hard pass in out_treatment1.bed and out_treatment2.bed as destinations
python3 ../scripts/wellington_bootstrap.py $@

#unlike the previous analysis, now everything is handily placed in the root

#do two visualisations, for each possible file
mkdir treatment1_output_visualisation
#get the wiggle tracks
python3 ../scripts/dnase_wig_tracks.py $3 $1 treatment1_output_visualisation/fw_cuts.wig treatment1_output_visualisation/rv_cuts.wig
if [ -s out_treatment1.bed ]
	then
		#get the average profile thing
		python3 ../scripts/dnase_average_profile.py out_treatment1.bed $1 treatment1_output_visualisation/average_footprint.png
		#get the heatmap
		python3 ../scripts/dnase_to_javatreeview.py out_treatment1.bed $1 treatment1_output_visualisation/javatreeview_heatmap_ready.csv
fi
mkdir treatment2_output_visualisation
#get the wiggle tracks
python3 ../scripts/dnase_wig_tracks.py $3 $2 treatment2_output_visualisation/fw_cuts.wig treatment2_output_visualisation/rv_cuts.wig
if [ -s out_treatment2.bed ]
	then
		#get the average profile thing
		python3 ../scripts/dnase_average_profile.py out_treatment1.bed $2 treatment2_output_visualisation/average_footprint.png
		#get the heatmap
		python3 ../scripts/dnase_to_javatreeview.py out_treatment1.bed $2 treatment2_output_visualisation/javatreeview_heatmap_ready.csv
fi

#get the wiggle tracks
python3 ../scripts/dnase_wig_tracks.py $1 $2 output_visualisation/fw_cuts.wig output_visualisation/rv_cuts.wig