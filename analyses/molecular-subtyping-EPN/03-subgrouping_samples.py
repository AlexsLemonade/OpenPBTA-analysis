# This script takes in the final EPN table with all the columns (in results folder  results/EPN_all.tsv), 
#	sets some rules and thresholds so samples can be  subgrouped
#	and adds the sub group name in the last "subgroup" column


## Importing modules to be used
import pandas as pd
import numpy as np 
import sys
import argparse
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--final_table', required = True,
                    help = 'path to the notebook with all EPN columns')
parser.add_argument('--subgroup_table', required = True,
		    help = 'output table  with subgroup')
parser.add_argument('--temp_folder_name', required = True,
                    help = 'outout table  with subgroup')
args = parser.parse_args()

# REading in input table  
EPN_final = pd.read_csv(args.final_table, sep="\t")
# Adding a column with empty strings so subgroup name can be added
EPN_final["subgroup"] = ""

# File to write all the data to
outfile = open(args.subgroup_table, "w")

# This function takes ina  row, value in the column from EPN_final table, 
#	subgroupname that  should be returned if threshold is met, 
#	value of the  threshold that  should  be higher in the row
# It returns the subgroupname(which will be added to the final column)
#	Before returning the subgroup name, it also checks if there is a subgroupname and
#	appends to it  with a "," separator  
def subgroup_func(row, columnname, subgroupname, value):
    if row[columnname] > value:
        if row["subgroup"] == '':
            return(subgroupname)
        elif subgroupname in row["subgroup"].split(","):
            return(row["subgroup"])
        else:
            return(row["subgroup"]+","+subgroupname)
    else:
        return(row["subgroup"])

## There is no set thresholds for gene expression outliers, so creating histograms for each of them. 
## If gene expression Z-score is very high(outliers), that gene might be overexpressed in that sample and 
## helps the sample to be subgrouped accordingly?? 
    
    
#### Looking for  ST_EPN_RELA sub-group samples
# PTEN--TAS2R1 fusion column is all zeros 
EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="C11orf95--RELA", subgroupname="ST_EPN_RELA", value=0)
EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="LTBP3--RELA", subgroupname="ST_EPN_RELA", value=0)
EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="9p_loss", subgroupname="ST_EPN_RELA", value=0)
EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="9q_loss", subgroupname="ST_EPN_RELA", value=0)
# Cannot find clear outliers from the histogram plots for RELA Z-scores 
# One Z-score above 3 
fig = plt.figure(figsize=(20,10))
L1CAM_fig = plt.figure()
plt.hist(list(EPN_final["L1CAM_expr_zscore"])) 
plt.suptitle('L1CAM expression')
fig.savefig(args.temp_folder_name+"/L1CAM_hist.png", dpi=L1CAM_fig.dpi)
EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="L1CAM_expr_zscore", subgroupname="ST_EPN_RELA", value=3)


EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="C11orf95--YAP1", subgroupname="ST_EPN_YAP1", value=0)
EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="YAP1--MAMLD1", subgroupname="ST_EPN_YAP1", value=0)
EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="C11orf95--MAML2", subgroupname="ST_EPN_YAP1", value=0)
EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="YAP1--FAM118B", subgroupname="ST_EPN_YAP1", value=0)
EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="11q_loss", subgroupname="ST_EPN_YAP1", value=0)
EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="11q_gain", subgroupname="ST_EPN_YAP1", value=0)
# Two sample higher than 1.7, look like outliers based on histogram
ARL4D_fig = plt.figure()
plt.hist(list(EPN_final["ARL4D_expr_zscore"])) 
plt.suptitle('ARL4D expression')
fig.savefig(args.temp_folder_name+"/ARL4D_hist.png", dpi=ARL4D_fig.dpi)
EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="ARL4D_expr_zscore", subgroupname="ST_EPN_YAP1", value=1.7)
# One sample higher than 1.0, another sample at 0.28 also look like an outlier but a cutoff value of 0.25 maybe too small?? 
CLDN1_fig = plt.figure()
plt.hist(list(EPN_final["CLDN1_expr_zscore"])) 
plt.suptitle('CLDN1 expression')
fig.savefig(args.temp_folder_name+"/CLDN1_hist.png", dpi=CLDN1_fig.dpi)
EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="CLDN1_expr_zscore", subgroupname="ST_EPN_YAP1", value=1)



EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="1q_loss", subgroupname="PT_EPN_A", value=0)
# Anything above 2 Z-score looks like an outlier
TKTL1_fig = plt.figure()
plt.hist(list(EPN_final["TKTL1_expr_zscore"])) 
plt.suptitle('TKTL1 expression')
fig.savefig(args.temp_folder_name+"/TKTL1_hist.png", dpi=TKTL1_fig.dpi)
EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="TKTL1_expr_zscore", subgroupname="PT_EPN_A", value=2)
# Any score above 2.0? Two samples at 2.03 and 2.07. Not includng this one in the PF_EPN_A
CXorf67_fig = plt.figure()
plt.hist(list(EPN_final["CXorf67_expr_zscore"])) 
plt.suptitle('CXorf67 expression')
fig.savefig(args.temp_folder_name+"/CXorf67_hist.png", dpi=CXorf67_fig.dpi)
EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="CXorf67_expr_zscore", subgroupname="PT_EPN_A", value=2)



EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="6q_loss", subgroupname="PT_EPN_B", value=0)
EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="6p_loss", subgroupname="PT_EPN_B", value=0)
# Two sample with Z-scores above 3
IFT46_fig = plt.figure()
plt.hist(list(EPN_final["IFT46_expr_zscore"])) 
plt.suptitle('IFT46 expression')
fig.savefig(args.temp_folder_name+"/IFT46_hist.png", dpi=IFT46_fig.dpi)
EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="IFT46_expr_zscore", subgroupname="PT_EPN_B", value=3)
# Look at GPB1 expression Z-scores  
GPBP1_fig = plt.figure()
plt.hist(list(EPN_final["GPBP1_expr_zscore"])) 
plt.suptitle('GPBP1 expression')
fig.savefig(args.temp_folder_name+"/GPBP1_hist.png", dpi=GPBP1_fig.dpi)
fig.savefig('/Users/kogantit/Documents/OpenPBTA/OpenPBTA-analysis/analyses/molecular-subtyping-EPN/GPBP1_hist.png', dpi=GPBP1_fig.dpi)
EPN_final["subgroup"] = EPN_final.apply(subgroup_func, axis=1, columnname="GPBP1_expr_zscore", subgroupname="PT_EPN_B", value=2)


#EPN_final.to_csv("results/EPN_all_data_withsubgroup.tsv", index=None)
EPN_final.to_csv(outfile, sep="\t", header=True, index=False)
outfile.close()


