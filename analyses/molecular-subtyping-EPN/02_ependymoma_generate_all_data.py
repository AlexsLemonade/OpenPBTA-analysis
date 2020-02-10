#!/usr/bin/env python
# Author Teja Koganti (D3B)

import sys

import argparse
import pandas as pd
import numpy as np 
import zipfile
from scipy import stats

parser = argparse.ArgumentParser()
parser.add_argument('--notebook', required = True,
                    help = 'path to the notebook file')
parser.add_argument('--gistic', required = True,
                    help = "gistic zip file")
parser.add_argument('--subfile-gistic', required = True,
                    help = "subfile of the zip folder that contains broad values by arm")
parser.add_argument('--gsva', required = True,
                    help = "gsva scores file")
parser.add_argument('--expression', required = True,
                    help = "expression subset file")
parser.add_argument('--fusion', required = True,
                    help = "fusion results file")
parser.add_argument('--breakpoints', required = True,
                    help = "breakpoint summary file")
parser.add_argument('-o', '--outfile', required = True,
                    help = "output file")
args = parser.parse_args()

# File to write all the data in 
outfile = open(args.outfile, "w")

# Reading GISTIC broad_values file for CNA 
zip=zipfile.ZipFile(args.gistic)
CNA=pd.read_csv(zip.open(args.subfile_gistic), sep="\t")
CNA = CNA.set_index("Chromosome Arm")


# Reading in gene set enrichment analyses file for GSEA scores for NKKB pathway
gsva = pd.read_csv(args.gsva, sep="\t")
gsva_NFKB = gsva.loc[gsva['hallmark_name'] == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"]
gsva_NFKB = gsva_NFKB.set_index("Kids_First_Biospecimen_ID")

# Reading subset gene expression file
fpkm_df = pd.read_csv(args.expression, sep = "\t")
fpkm_df = fpkm_df.set_index("GENE")
zscore_fpkm_df = fpkm_df.apply(stats.zscore)

# Reading fusion summary file
fusion = pd.read_csv(args.fusion, sep="\t")
fusion = fusion.set_index("Kids_First_Biospecimen_ID")

# Reading chromosomal instability file for breakpoint density 
breakpoint_density = pd.read_csv(args.breakpoints, sep="\t")
breakpoint_density = breakpoint_density.set_index("samples")

# Reading the input in a  dataframe 
EPN_notebook = pd.read_csv(args.notebook, sep="\t")

# This function takes in CNA dataframe along  with chromosomal arm
# and if should be a loss/gain in the results and returns a boolean value  accordingly
def DNA_samples_fill_df(row, CNA, arm, loss_gain):
    if row["Kids_First_Biospecimen_ID_DNA"] is np.nan:
        return(0)
    CNA_value = CNA.loc[arm, row["Kids_First_Biospecimen_ID_DNA"]]
    if CNA_value<0 and loss_gain=="loss":
       return(1)
    elif loss_gain=="gain" and CNA_value>0:
       return(1)
    else:
       return(0)

# Function to generate Z-scores column for every gene 
def fill_df_with_fpkm_zscores(df, fpkmdf, gene_name):
        zscore_list = stats.zscore(np.array(df.apply(lambda x: fpkmdf.loc[gene_name, x["Kids_First_Biospecimen_ID_RNA"]], axis=1)))
        column_name = gene_name + "_expr_zscore"
        df[column_name] = pd.Series(zscore_list)
        return(df)


# This  function takes a dataframe whose values need to be  used for final 
# EPN_notebook and based on columnname, the corresponding value is returned  
def  fill_df_with_rnaresults(row, RNA_df, columnname):
	fusion_value = RNA_df.loc[row["Kids_First_Biospecimen_ID_RNA"],  columnname]
	return(fusion_value)

# Filling up  dataframe with broad CNA values from GISTIC
EPN_notebook["1q_loss"] = EPN_notebook.apply(DNA_samples_fill_df, axis = 1, args = (CNA, "1q", "loss"))
EPN_notebook["9p_loss"] = EPN_notebook.apply(DNA_samples_fill_df, axis = 1, args = (CNA, "9p", "loss"))
EPN_notebook["9q_loss"] = EPN_notebook.apply(DNA_samples_fill_df, axis = 1, args = (CNA, "9q", "loss"))
EPN_notebook["6p_loss"] = EPN_notebook.apply(DNA_samples_fill_df, axis = 1, args = (CNA, "6p", "loss"))
EPN_notebook["6q_loss"] = EPN_notebook.apply(DNA_samples_fill_df, axis = 1, args = (CNA, "6q", "loss"))
EPN_notebook["11q_loss"] = EPN_notebook.apply(DNA_samples_fill_df, axis = 1, args = (CNA, "11q", "loss"))
EPN_notebook["11q_gain"] = EPN_notebook.apply(DNA_samples_fill_df, axis = 1, args = (CNA, "11q", "gain"))

# Adding GSEA score to the  dataframe 
#EPN_notebook["NFKB_pathway_GSEAscore"] = EPN_notebook.apply(lambda x: gsva_NFKB.loc[x["Kids_First_Biospecimen_ID_RNA"], "gsea_score"], axis =1)
EPN_notebook["NFKB_pathway_GSEAscore"] = EPN_notebook.apply(fill_df_with_rnaresults, axis=1, args =  (gsva_NFKB, "gsea_score"))

fusions_list = ["C11orf95--RELA", "LTBP3--RELA", "PTEN--TAS2R1",  "C11orf95--YAP1", "YAP1--MAMLD1", "YAP1--FAM118B", "C11orf95--MAML2"]
for every_fusion in fusions_list:
	EPN_notebook[every_fusion] = EPN_notebook.apply(fill_df_with_rnaresults,  axis=1, args = (fusion, every_fusion))


# Adding breakpoints density for chromosomal instability  to the dataframe
EPN_notebook["breaks_density-chromosomal_instability"] = EPN_notebook.apply(lambda x: breakpoint_density.loc[x["Kids_First_Biospecimen_ID_DNA"], "breaks_density"]
        if x["Kids_First_Biospecimen_ID_DNA"] is not np.nan  else "NA", axis=1)

# Adding Z-scores to dataframe
EPN_notebook = fill_df_with_fpkm_zscores(EPN_notebook, fpkm_df, "RELA")
EPN_notebook = fill_df_with_fpkm_zscores(EPN_notebook, fpkm_df, "L1CAM")
EPN_notebook = fill_df_with_fpkm_zscores(EPN_notebook, fpkm_df, "ARL4D")
EPN_notebook = fill_df_with_fpkm_zscores(EPN_notebook, fpkm_df, "CLDN1")
EPN_notebook = fill_df_with_fpkm_zscores(EPN_notebook, fpkm_df, "CXorf67")
EPN_notebook = fill_df_with_fpkm_zscores(EPN_notebook, fpkm_df, "TKTL1")
EPN_notebook = fill_df_with_fpkm_zscores(EPN_notebook, fpkm_df, "GPBP1")
EPN_notebook = fill_df_with_fpkm_zscores(EPN_notebook, fpkm_df, "IFT46")

# Replacing all Nan values with NA so they are not empty when writing to a file
EPN_notebook =EPN_notebook.replace(np.nan, 'NA', regex=True)
# Writing dataframe to output file
EPN_notebook.to_csv(outfile, sep="\t", header=True, index=False)
outfile.close()

