# Chante Bethell for CCDL 2019 and Jo Lynne Rokita
# Code adapted from the PPTC PDX Oncoprint Generation repository here:
# https://github.com/marislab/create-pptc-pdx-oncoprints/tree/master/R
#
# # #### USAGE
# This script is intended to be sourced in the 
# 'analyses/oncoprint-landscape/01-plot-oncoprint.R' script as follows:
# 
# source(file.path("util", "oncoplot-palette.R"))


# Define a color vector for plots 
color_palette <- c("Missense_Mutation" = "#35978f", 
                   "Nonsense_Mutation" = "#000000",
                   "Frame_Shift_Del" = "#56B4E9", 
                   "Frame_Shift_Ins" = "#FFBBFF", 
                   "Splice_Site" = "#F0E442",
                   "Translation_Start_Site" = "#191970",
                   "Nonstop_Mutation" = "#545454",
                   "In_Frame_Del" = "#CAE1FF",
                   "In_Frame_Ins" = "#FFE4E1",
                   "Stop_Codon_Ins" = "#CC79A7",
                   "Start_Codon_Del" = "#56B4E9",
                   "Fusion" = "#7B68EE",
                   "Multi_Hit" = "#f46d43",
                   "Hom_Deletion" = "#313695",
                   "Hem_Deletion" = "#abd9e9",
                   "Amplification" = "#c51b7d",
                   "Loss" = "#0072B2",
                   "Gain" = "#D55E00",
                   "High_Level_Gain" = "#FF0000",
                   "Multi_Hit_Fusion" = "#CD96CD")