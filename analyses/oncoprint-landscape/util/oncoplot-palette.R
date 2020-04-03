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
                   "Multi_Hit" = "#00F021",
                   "Hom_Deletion" = "#313695",
                   "Hem_Deletion" = "#abd9e9",
                   "amplification" = "#c51b7d",
                   "loss" = "#0072B2",
                   "gain" = "#D55E00",
                   "High_Level_Gain" = "#FF0000",
                   "Multi_Hit_Fusion" = "#CD96CD",
                   "Adenoma" = "#f23d3d",
                   "ATRT" =	"#731d1d",
                   "Central neurocytoma" = "#b38686",
                   "Chondrosarcoma" =	"#cc5c33",
                   "Chordoma" =	"#331c0d",
                   "Choroid plexus tumor" =	"#ffb380",
                   "CNS EFT-CIC" = "#b25f00",
                   "CNS lymphoma"	= "#f2d6b6",
                   "CNS neuroblastoma" = "#736556",
                   "CNS Rhabdomyosarcoma"	= "#ffaa00",
                   "CNS sarcoma" =	"#4c3d00",
                   "Craniopharyngioma" =	"#e2f200",
                   "DNET" =	"#919926",
                   "Dysplasia" = "#d6f2b6",
                   "Embryonal Tumor" = "#304d26",
                   "Ependymoma" = "#00f241",
                   "ETMR" =	"#009929",
                   "Ganglioglioma" = "#698c7c",
                   "Germinoma" = "#39e6c3",
                   "Glial-neuronal tumor NOS" = "#005359",
                   "Gliosis" =	"#263233",
                   "Hemangioblastoma" =	"#00c2f2",
                   "Hemangioma"	= "#40a6ff",
                   "HGAT" =	"#406280",
                   "Langerhans Cell histiocytosis" = "#0044ff",
                   "LGAT" =	"#00144d",
                   "LGMT" = "#acbbe6",
                   "Medulloblastoma" = "#7373e6",
                   "Meningioma" =	"#3d0099",
                   "MPNST" = "#c200f2",
                   "Neurofibroma" =	"#917399",
                   "none" =	"#f1f1f1",
                   "Oligodendroglioma" = "#f279da",
                   "Other" = "#cc0052",
                   "Pineoblastoma" = "#994d6b",
                   "Schwannoma" =	"#4d2636",
                   "Teratoma" =	"#ffbfd9"
)
