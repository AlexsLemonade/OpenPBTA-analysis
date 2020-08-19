##  This script takes in TMB scores outfile and calculates cumulative distribution
#   function plot.
# INputs  include the tmb  scores text file, outfile name and min number of samples to be plotted
#                   (no need to use any extensions for out file, just the name)
#   Somtimes MAF can have very few samples under each histology and those histologies can be filtered
#       out of the plot using the min samples feature

# The TMB  text  should have the following column headers -
#   Samplename experimental_strategy disease count bedlength  TMB
# Use results from 01_calculate_tmb_targetflexible.py script

####### Author: Teja Koganti ##################
# python3 02_cumulative_freq_TMBplot.py  \
#   -i ../output/pbta-snv-consensus.tmb.txt \
#   -o ../output/pbta-snv-consensus \
#   -s 10

import argparse
import subprocess
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--tmb_scores", required=True, help="file with TMB scores")
parser.add_argument(
    "-o", "--outfilename", required=True, help="Name of the out plot, no extension"
)
parser.add_argument(
    "-s",
    "--minsamplestoplot",
    required=True,
    help="Minimum samples from each histology/disease to plot",
)
args = parser.parse_args()

######### Preparing inputs files and disease list ###########
# Loading TMB scres text file into a  dataframe
tmbfile = pd.read_csv(args.tmb_scores, sep="\t")

# Choosing disease list to plot
#tmbfile["new_TMB"] = ( tmbfile["count"] + 1 * 1000000) / tmbfile["bedlength"]
#tmbfile["new_TMB"] = tmbfile["count"] + 1 / tmbfile["bedlength"]


tmbfile["new_TMB"] = tmbfile.apply(lambda x:
    (x["count"] + 1 * 1000000) / x["bedlength"] if x["count"]==0
    else x["TMB"],axis=1)
#( tmbfile["count"] + 1 * 1000000) / tmbfile["bedlength"]
#value_to_add = np.min(tmbfile[tmbfile["TMB"] != 0]["TMB"])*0.8
#tmbfile["new_TMB"] = tmbfile["TMB"].replace({0: value_to_add})
#tmbfile = tmbfile.sort_values(by=['new_TMB'])



# Counting number of lines under each disease and saving it as disease_number
# If lines per disease is greater than minsamplestoplot, then add disease to top_diseases

# Count how many samples in each disease
disease_number = tmbfile.groupby("disease").count()["new_TMB"]
# List only diseases that have a higher number of samples than minsamples input
top_diseases = disease_number[disease_number > int(args.minsamplestoplot)].index
# Sort disease list based on median value. Output plot will use this sorting
sorted_diseases = tmbfile.groupby("disease").median().sort_values("new_TMB")


sorted_diseases_final = []
for dis in sorted_diseases.index:
    if dis in top_diseases:
        sorted_diseases_final.append(dis)

##############################################################


############### Creating TMB plots #############################
# Creating TMB cumulative dist plot
perdisease = tmbfile.groupby("disease")
plt.figure(figsize=(20, 10))
tick_locs = []
for i, disease in enumerate(
    sorted_diseases_final
):  # i is index number, disease is disaese name
    disease_df = perdisease.get_group(disease)  # Getting disease DF
    TMB = disease_df["new_TMB"].sort_values()  # TMB will ne used for y-axis
    x = (
        np.arange(1, len(TMB) + 1) / len(TMB)
    ) + i  # this generate an array in every loop that clusters
    # all disease x-axis together and adds so the next disease
    # x-axis numbers are spaced separately from previous disease
    # print(np.median(x))
    tick_locs.append(
        np.median(x)
    )  # Takes the median of x-axis array and prints disease name there
    # print((i+max(x))/2)
    plt.plot(x, TMB, linestyle="none", marker="o")
    plt.hlines(np.median(TMB), np.percentile(x, 30), np.percentile(x, 70))
    # Takes y-axis median to know where to print median line
    # 30 and 70 show how long the median line should go on x-axis
    plt.yscale("log")
plt.xticks(tick_locs, rotation=90, fontsize=15)
ax = plt.gca()
ax.set_xticklabels(sorted_diseases_final)
ax.set_ylabel("Mutations per Mb", fontsize=20)
ax.set_xlabel("Disease type", fontsize=20)
plt.tight_layout()
plt.savefig(args.outfilename + ".tmb.png")
###############################################################
