#!/usr/bin/env python
# Author Teja Koganti (D3B)

#import sys
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
parser.add_argument('--subfile-gistic-broad', required = True,
                    help = "subfile of the zip folder that contains broad values by arm")
parser.add_argument('--gsva', required = True,
                    help = "gsva scores file")
parser.add_argument('--expression', required = True,
                    help = "expression subset file")
parser.add_argument('--fusion', required = True,
                    help = "fusion results file")
parser.add_argument('--breakpoints-cnv', required = True,
                    help = "breaks density  CNV summary file")
parser.add_argument('--breakpoints-sv', required = True,
                    help = "breaks density  SV summary file")
parser.add_argument('-o', '--outfile', required = True,
                    help = "output file")
parser.add_argument('--focal-gene-cn', required = True,
                    help = "focal-cn-file-preparation based on gene")
parser.add_argument('--subfile-gistic-focalbygene', required = True,
                    help = "subfile of the GISTIC zip folder that contains focal data fro CDKN2A")
args = parser.parse_args()

# File to write all the data in
outfile = open(args.outfile, "w")

# Reading GISTIC broad_values and focal_by_genefile for CNA
gisticzip = zipfile.ZipFile(args.gistic)
broad_CNA = pd.read_csv(gisticzip.open(args.subfile_gistic_broad), sep="\t")
broad_CNA = broad_CNA.set_index("Chromosome Arm")

gistic_focalCN = pd.read_csv(gisticzip.open(args.subfile_gistic_focalbygene), sep="\t")
gistic_focalCN = gistic_focalCN.set_index("Gene Symbol")
#Select only the column we need and transform to make consistent: sample names in index
gistic_focalCN = gistic_focalCN.loc[["CDKN2A",]].T


# Reading in gene set enrichment analyses file for GSEA scores for NKKB pathway
gsva = pd.read_csv(args.gsva, sep="\t")
gsva_NFKB = gsva.loc[gsva['hallmark_name'] == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"]
gsva_NFKB = gsva_NFKB.set_index("Kids_First_Biospecimen_ID")

# Reading subset gene expression file
fpkm_df = pd.read_csv(args.expression, sep = "\t")
fpkm_df = fpkm_df.set_index("GENE")

# Reading fusion summary file
fusion_df = pd.read_csv(args.fusion, sep="\t")
fusion_df = fusion_df.set_index("Kids_First_Biospecimen_ID")

# Reading chromosomal instability file for breakpoint density for CNV
breakpoint_density_cnv = pd.read_csv(args.breakpoints_cnv, sep="\t")
breakpoint_density_cnv = breakpoint_density_cnv.set_index("samples")

# Reading chromosomal instability  file for breakpoint for SV
breakpoint_density_sv  = pd.read_csv(args.breakpoints_sv, sep="\t")
breakpoint_density_sv  = breakpoint_density_sv.set_index("samples")

# Reading consensus focal CN results  from analyses
focal_cn_gene = pd.read_csv(args.focal_gene_cn, sep="\t")
focal_cn_gene_CDKN2A = focal_cn_gene.loc[focal_cn_gene['gene_symbol'] == "CDKN2A"]
focal_cn_gene_CDKN2A = focal_cn_gene_CDKN2A.set_index("biospecimen_id")

# Reading the input in a  dataframe
EPN_notebook = pd.read_csv(args.notebook, sep="\t",index_col=False)



# This  function takes a dataframe whose values need to be  used for final
# EPN_notebook and based on row_name, the corresponding value is returned
# If a sample is not in included_samples, it's values are set to NA
# (If included samples is blank, this is ignored)
def fill_df(sample, DNA_df, col_name, included_samples = None, default = 0):
    if sample is np.nan: # no results for this sample
        return("NA")
    elif included_samples and sample not in included_samples:
        return("NA")
    elif sample not in DNA_df.index.to_list():
        # sample is not present in the data set, but should get a value
        return(default)
    else:
        value = DNA_df.loc[sample, col_name]
        return(value)

# This function takes in CNA dataframe along  with chromosomal arm
# and if should be a loss/gain in the results and returns a boolean value accordingly
def broad_CNA_fill_df(row, CNA, arm, loss_gain):
    if row["Kids_First_Biospecimen_ID_DNA"] is np.nan: # no DNA results for this sample
        return("NA")
    CNA_value = CNA.loc[arm, row["Kids_First_Biospecimen_ID_DNA"]]
    if CNA_value is np.nan:
        return ("NA")
    elif CNA_value < 0 and loss_gain=="loss":
       return("1")
    elif loss_gain=="gain" and CNA_value>0:
       return("1")
    else:
       return("0")


# Function to generate Z-scores column for every gene
def fill_df_with_fpkm_zscores(df, fpkmdf, gene_name):
    zscore_list = stats.zscore(np.array(df.loc[df["Kids_First_Biospecimen_ID_RNA"].notna(),:].apply(lambda x: fpkmdf.loc[gene_name, x["Kids_First_Biospecimen_ID_RNA"]], axis=1)))
    column_name = gene_name + "_expr_zscore"
    # add z-score array to df_rna column_name
    df.loc[df["Kids_First_Biospecimen_ID_RNA"].notna(),column_name] = pd.Series(zscore_list,index=df[df["Kids_First_Biospecimen_ID_RNA"].notna()].index.array)
    # add NA to expression zscore columns for dna only samples
    df.loc[df["Kids_First_Biospecimen_ID_RNA"].isna(),column_name] = np.nan

    return(df)



# Get the list of DNA samples that made it through the pipeline
# All that are present in the GISTIC data passed consensus filtering.
# Others may be included or excluded in other CN data sets, but should be set to NA
cn_called_samples = broad_CNA.columns.to_list()

# Filling up  dataframe with broad CNA values from GISTIC
EPN_notebook["1q_loss"]  = EPN_notebook.apply(broad_CNA_fill_df, axis = 1, args = (broad_CNA, "1q",  "loss"))
EPN_notebook["1q_gain"]  = EPN_notebook.apply(broad_CNA_fill_df, axis = 1, args = (broad_CNA, "1q",  "gain"))
EPN_notebook["9p_loss"]  = EPN_notebook.apply(broad_CNA_fill_df, axis = 1, args = (broad_CNA, "9p",  "loss"))
EPN_notebook["9q_loss"]  = EPN_notebook.apply(broad_CNA_fill_df, axis = 1, args = (broad_CNA, "9q",  "loss"))
EPN_notebook["6p_loss"]  = EPN_notebook.apply(broad_CNA_fill_df, axis = 1, args = (broad_CNA, "6p",  "loss"))
EPN_notebook["6q_loss"]  = EPN_notebook.apply(broad_CNA_fill_df, axis = 1, args = (broad_CNA, "6q",  "loss"))
EPN_notebook["11q_loss"] = EPN_notebook.apply(broad_CNA_fill_df, axis = 1, args = (broad_CNA, "11q", "loss"))
EPN_notebook["11q_gain"] = EPN_notebook.apply(broad_CNA_fill_df, axis = 1, args = (broad_CNA, "11q", "gain"))

# Adding GSEA score to the  dataframe
#EPN_notebook["NFKB_pathway_GSEAscore"] = EPN_notebook.apply(lambda x: gsva_NFKB.loc[x["Kids_First_Biospecimen_ID_RNA"], "gsea_score"], axis =1)
EPN_notebook["NFKB_pathway_GSEAscore"] = EPN_notebook["Kids_First_Biospecimen_ID_RNA"].apply(
    fill_df,
    args = (gsva_NFKB, "gsea_score")
    )

# Adding all fusions to a list so they can use the function
# fill_df to fill appropriate fusion summary information under each fusion
fusions_list = ["C11orf95--RELA", "LTBP3--RELA", "PTEN--TAS2R1",  "C11orf95--YAP1", "YAP1--MAMLD1", "YAP1--FAM118B", "C11orf95--MAML2", "YAP1--MAML2"]
for fusion in fusions_list:
	EPN_notebook[fusion] = EPN_notebook["Kids_First_Biospecimen_ID_RNA"].apply(
        fill_df,
        args = (fusion_df, fusion)
        )


# Adding breakpoints density for chromosomal instability  to the dataframe
EPN_notebook["breaks_density-chromosomal_instability_CNV"] = EPN_notebook["Kids_First_Biospecimen_ID_DNA"].apply(
    fill_df,
    args = (breakpoint_density_cnv, "breaks_density", cn_called_samples)
    )
EPN_notebook["breaks_density-chromosomal_instability_SV"]  = EPN_notebook["Kids_First_Biospecimen_ID_DNA"].apply(
    fill_df,
    args = (breakpoint_density_sv,  "breaks_density", cn_called_samples)
    )

# Adding focal CN from GISTIC files for CDKN2A
EPN_notebook["GISTIC_focal_CN_CDKN2A"] = EPN_notebook["Kids_First_Biospecimen_ID_DNA"].apply(
    fill_df,
    args = (gistic_focalCN, "CDKN2A", cn_called_samples)
    )

# Adding focal CN from CDKN2A from CNV consensus files in analyses
# Using status column from consensus_seg_annotated_cn_autosomes.tsv.gz file
EPN_notebook["consensus_focal_CN_CDKN2"] = EPN_notebook["Kids_First_Biospecimen_ID_DNA"].apply(
    fill_df,
    args = (focal_cn_gene_CDKN2A, "status", cn_called_samples)
    )

# Adding Z-scores to dataframe
expression_cols = ["RELA", "L1CAM", "ARL4D", "CLDN1", "CXorf67", "TKTL1", "GPBP1", "IFT46"]

for gene in expression_cols:
    # rna matched dna samples
    EPN_notebook = fill_df_with_fpkm_zscores(EPN_notebook, fpkm_df, gene)
        
# Replacing all Nan values with NA so they are not empty when writing to a file
EPN_notebook = EPN_notebook.replace(np.nan, 'NA', regex=True)
# Sort
EPN_notebook = EPN_notebook.sort_values(by = ["Kids_First_Participant_ID", "sample_id"])

# Writing dataframe to output file
EPN_notebook.to_csv(outfile, sep="\t", header=True, index=False)
outfile.close()

