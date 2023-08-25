#!/bin/sh

# print start time
echo "start: " `date`


########################################################################
#                         Info on script etc                           #
########################################################################

# Downloaded nsegata-lefse-9adc3a62460e folder from https://bitbucket.org/nsegata/lefse/downloads/

# How I made the earlypython virtual environment:
# conda create --name earlypython python=2.7
# conda activate earlypython
# conda info --envs
# python --version
# conda install -c bioconda lefse #success!! It wasn't working on the current python version

# On the command line you can use [conda activate /Users/aliciahadingham/opt/miniconda3/envs/earlypython] to activate the early python virtual environment. 
# But in this sh script only [source /Users/aliciahadingham/anaconda3/etc/profile.d/conda.sh] [conda activate /Users/aliciahadingham/opt/miniconda3/envs/earlypython] works. 
# Help from: https://github.com/conda/conda/issues/7980



########################################################################
#                          Prep for LEfSe                              #
########################################################################

# Change wd to the one where we will run LEfSe in
cd /Users/aliciahadingham/Desktop/nsegata-lefse-9adc3a62460e

# Remove all previously generated files and plots, so that we know only files in these folders relate to just this run of LEfSe
rm data_for_lefse/*
rm output_RES/*
rm output_lda_plots/*
rm output_cladograms/*

# Change wd
cd /Users/aliciahadingham/OneDrive\ -\ King\'s\ College\ London/PhD/Projects/BV_vs_OTUs/data_BV/LEfSe

# Copying all txt files (files needed for LEfSe that were exported from R) in the wd across
cp -R *.txt /Users/aliciahadingham/Desktop/nsegata-lefse-9adc3a62460e/data_for_lefse

# Activate earlypython virtual environment
source /Users/aliciahadingham/anaconda3/etc/profile.d/conda.sh
conda activate /Users/aliciahadingham/opt/miniconda3/envs/earlypython

# Checking we are in the earlypython virtual environment
conda info --envs
# This should have an * by earlypython 

# Checking python version in this virtual environment
python --version
# Python 2.7.15

# Install necessary packages if not yet installed:
# conda install -c biobakery lefse
# conda install -c r rpy2

# Look at help on the python scripts:
# python format_input.py --help
# python run_lefse.py --help
# python plot_res.py --help
# python plot_cladogram.py --help

# Change wd to the one with the lefse python scripts
cd /Users/aliciahadingham/Desktop/nsegata-lefse-9adc3a62460e




########################################################################
#                   Formating input file for LEfSe                     #       
########################################################################

# TXT file path
TXTs_from_R="data_for_lefse/*.txt"

# Command to format the input file (file.txt) & generate a output file (file.in)
for i in $TXTs_from_R
do
python format_input.py $i ${i/.txt/.in} -c 1 -u 2 -o 1000000
done
# Arguments:
# -c [1..n_feats]       set which feature use as class (default 1)
# -s [1..n_feats]       set which feature use as subclass (default -1 meaning no subclass)
# -u [1..n_feats]       set which feature use as subject (default -1 meaning no subject)
# -o float              set the normalization value (default -1.0 meaning no normalization)      



########################################################################
#                           Running LEfSe                              #
########################################################################

# IN file path - effect of BV
lefse_input__BV_status="data_for_lefse/BV_status*.in"

# IN file path - effect of ethnicity
lefse_input__ethnicity="data_for_lefse/Ethnicity*.in"


# We are going to use a LDA cut off of 4 for looking at BV status and 4 for looking at ethnicity.

# Generate a file (file.res) of the LEfSe analysis results, based on the input file (file.in)
for i in $lefse_input__BV_status
do
tmp=${i/data_for_lefse/output_RES}
python run_lefse.py $i ${tmp/.in/.res} -l 4
done
# Arguments:
# -l float        set the threshold on the absolute value of the logarithmic LDA score (default 2.0)


# Generate a file (file.res) of the LEfSe analysis results, based on the input file (file.in)
for i in $lefse_input__ethnicity
do
tmp=${i/data_for_lefse/output_RES}
python run_lefse.py $i ${tmp/.in/.res} -l 4
done


########################################################################
#                      Plottting LEfSe results                         #
########################################################################

# RES file path  - BV vs normal
lefse_results__BV_vs_normal="output_RES/BV_status*BV_vs_normal*.res"

# RES file path  - BV vs int
lefse_results__BV_vs_int="output_RES/BV_status*BV_vs_int*.res"

# RES file path  - int vs normal
lefse_results__int_vs_normal="output_RES/BV_status*int_vs_normal*.res"

# RES file path  - effect of ethnicity
lefse_results__ethnicity="output_RES/Ethnicity*.res"


# We are going to run the lefse results looking at BV status over the python script where the colours are red, orange & green. The ethnicity results will be run over the python script where the colours will be blue and purple.

# BV_vs_normal
# To plot the results of the LEfSe analysis generated from the previous step, run the following command (first line). To visualize the results in a Cladogram, run the following command to generate the Cladogram figure (second line). 
for i in $lefse_results__BV_vs_normal
do
tmp1=${i/output_RES/output_lda_plots}
tmp2=${i/output_RES/output_cladograms}
python plot_res_red-green.py --format pdf $i ${tmp1/.res/.pdf}
python plot_cladogram_red-green.py --dpi 600 --format png $i ${tmp2/.res/_cladogram.png}
done

# BV_vs_int
# To plot the results of the LEfSe analysis generated from the previous step, run the following command (first line). To visualize the results in a Cladogram, run the following command to generate the Cladogram figure (second line). 
for i in $lefse_results__BV_vs_int
do
tmp1=${i/output_RES/output_lda_plots}
tmp2=${i/output_RES/output_cladograms}
python plot_res_red-orange.py --format pdf $i ${tmp1/.res/.pdf}
python plot_cladogram_red-orange.py --dpi 600 --format png $i ${tmp2/.res/_cladogram.png}
done

# int_vs_normal
# To plot the results of the LEfSe analysis generated from the previous step, run the following command (first line). To visualize the results in a Cladogram, run the following command to generate the Cladogram figure (second line). 
for i in $lefse_results__int_vs_normal
do
tmp1=${i/output_RES/output_lda_plots}
tmp2=${i/output_RES/output_cladograms}
python plot_res_orange-green.py --format pdf $i ${tmp1/.res/.pdf}
python plot_cladogram_orange-green.py --dpi 600 --format png $i ${tmp2/.res/_cladogram.png}
done



# Ethnicity
# To plot the results of the LEfSe analysis generated from the previous step, run the following command (first line). To visualize the results in a Cladogram, run the following command to generate the Cladogram figure (second line). 
for i in $lefse_results__ethnicity
do
tmp1=${i/output_RES/output_lda_plots}
tmp2=${i/output_RES/output_cladograms}
python plot_res_blue-magenta.py --format pdf $i ${tmp1/.res/.pdf}
python plot_cladogram_blue-magenta.py --dpi 600 --format png $i ${tmp2/.res/_cladogram.png}
done


########################################################################
#                               The End                                #
########################################################################

# print end time 
echo "end: " `date`


# How to run this script (change to the script's directory and then run using bash command):
# cd /Users/aliciahadingham/OneDrive\ -\ King\'s\ College\ London/PhD/Projects/BV_vs_OTUs/scripts
# bash LEfSe_BV_OTUs_with_loops_to_run_over_ALL_txts.sh
