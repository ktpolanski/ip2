#!/bin/bash

#pre-process peaks
cut -f -3 $3 > temp.bed
../bedtools/bin/bedops --range 50 --everything temp.bed > temp2.bed
../bedtools/bin/bedops --merge temp2.bed > processed_peaks.bed
rm temp.bed
rm temp2.bed

#run wellington.
#we need to run the indexing ourselves because yes
samtools index $1
samtools index $2
#hard pass in out_treatment1.bed and out_treatment2.bed as destinations
python ../scripts/wellington_bootstrap.py $1 $2 processed_peaks.bed out_treatment1.bed out_treatment2.bed ${@:4}

#unlike the previous analysis, now everything is handily placed in the root

#do two visualisations, for each possible file
mkdir treatment1_output_visualisation
#get the wiggle tracks
python ../scripts/dnase_wig_tracks.py processed_peaks.bed $1 treatment1_output_visualisation/fw_cuts.wig treatment1_output_visualisation/rv_cuts.wig
if [ -s out_treatment1.bed ]
	then
		#get the average profile thing
		python ../scripts/dnase_average_profile.py out_treatment1.bed $1 treatment1_output_visualisation/average_footprint.png
		#get the heatmap
		python ../scripts/dnase_to_javatreeview.py out_treatment1.bed $1 treatment1_output_visualisation/javatreeview_heatmap_ready.csv
fi
mkdir treatment2_output_visualisation
#get the wiggle tracks
python ../scripts/dnase_wig_tracks.py processed_peaks.bed $2 treatment2_output_visualisation/fw_cuts.wig treatment2_output_visualisation/rv_cuts.wig
if [ -s out_treatment2.bed ]
	then
		#get the average profile thing
		python ../scripts/dnase_average_profile.py out_treatment1.bed $2 treatment2_output_visualisation/average_footprint.png
		#get the heatmap
		python ../scripts/dnase_to_javatreeview.py out_treatment1.bed $2 treatment2_output_visualisation/javatreeview_heatmap_ready.csv
fi

#get the wiggle tracks
python ../scripts/dnase_wig_tracks.py $1 $2 output_visualisation/fw_cuts.wig output_visualisation/rv_cuts.wig