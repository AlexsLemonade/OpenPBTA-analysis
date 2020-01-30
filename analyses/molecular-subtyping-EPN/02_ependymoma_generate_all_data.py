#!/usr/bin/env python
# coding: utf-8
# Author Teja Koganti (D3B)


import pandas as pd
import numpy as np 
import zipfile
import  pyreadr
import statistics
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()

# File to write all the data in analyses/molecul-EPN/results folder 
outfile = open("analyses/molecular-subtyping-EPN/results/EPN_all_data.tsv", "w")

# Reading GISTIC broad_values file for CNA 
#CNA = pd.read_csv("/Users/kogantit/Documents/OpenPBTA_tickets/ependymoma_subtyping_ticket/2019-12-10-gistic-results-cnvkit/broad_values_by_arm.txt", sep="\t")
#CNA = CNA.set_index("Chromosome Arm")
zip=zipfile.ZipFile('data/pbta-cnv-cnvkit-gistic.zip')
CNA=pd.read_csv(zip.open('2019-12-10-gistic-results-cnvkit/broad_values_by_arm.txt'), sep="\t")
CNA = CNA.set_index("Chromosome Arm")

count=0

# Reading in gene set enrichment analyses file for GSEA scores for NKKB pathway
gsva = pd.read_csv("analyses/gene-set-enrichment-analysis/results/gsva_scores_stranded.tsv", sep="\t")
gsva_NFKB = gsva.loc[gsva['hallmark_name'] == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"]
gsva_NFKB = gsva_NFKB.set_index("Kids_First_Biospecimen_ID")

# Reading collapsed gene expression file, none in polya file, so ignoring that file 
readRDS = robjects.r["readRDS"]
fpkm_df = readRDS("analyses/molecular-subtyping-EPN/epn-subset/epn_pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds")

#Reading fusion summary file
fusion = pd.read_csv("analyses/fusion-summary/results/fusion_summary_ependymoma_foi.tsv", sep="\t")
fusion = fusion.set_index("Kids_First_Biospecimen_ID")

#Reading chromosomal instability file for breakpoint density 
breakpoint_density = pd.read_csv("analyses/chromosomal-instability/breakpoint-data/union_of_breaks_densities.tsv", sep="\t")
breakpoint_density = breakpoint_density.set_index("samples")

# creating empty lists to get mean and standard deviation of gene FPKM's from pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds file
#BSID_list = []
RELA_fpkm = []
L1CAM_fpkm = []
ARL4D_fpkm = []
CLDN1_fpkm = []
CXorf67_fpkm = []
TKTL1_fpkm = []
GPBP1_fpkm = []
IFT46_fpkm = []


# Opening the FPKM file separately to fill up the mean and stdev lists 
with open("analyses/molecular-subtyping-EPN/results/EPN_molecular_subtype.tsv", "r") as notebook:
    for line in notebook:
        line = line.split()
        if not line[0].startswith("Kids"):
            RELA_fpkm.append(fpkm_df.get_value("RELA", line[3]))
            L1CAM_fpkm.append(fpkm_df.get_value("L1CAM", line[3]))
            ARL4D_fpkm.append(fpkm_df.get_value("ARL4D", line[3]))
            CLDN1_fpkm.append(fpkm_df.get_value("CLDN1", line[3]))
            CXorf67_fpkm.append(fpkm_df.get_value("CXorf67", line[3]))
            TKTL1_fpkm.append(fpkm_df.get_value("TKTL1", line[3]))
            GPBP1_fpkm.append(fpkm_df.get_value("GPBP1", line[3]))
            IFT46_fpkm.append(fpkm_df.get_value("IFT46", line[3]))                  
                                

# Generating mean and stdev scores from the FPKM list
RELA_mean = np.mean(RELA_fpkm)
RELA_stdev = np.std(RELA_fpkm)
L1CAM_mean = np.mean(L1CAM_fpkm)
L1CAM_stdev = np.std(L1CAM_fpkm)
ARL4D_mean = np.mean(ARL4D_fpkm)
ARL4D_stdev = np.std(ARL4D_fpkm)
CLDN1_mean = np.mean(CLDN1_fpkm)
CLDN1_stdev = np.std(CLDN1_fpkm)
CXorf67_mean = np.mean(CXorf67_fpkm)
CXorf67_stdev = np.std(CXorf67_fpkm)
TKTL1_mean = np.mean(TKTL1_fpkm)
TKTL1_stdev = np.std(TKTL1_fpkm)
GPBP1_mean = np.mean(GPBP1_fpkm)
GPBP1_stdev = np.std(GPBP1_fpkm)
IFT46_mean = np.mean(IFT46_fpkm)
IFT46_stdev = np.std(IFT46_fpkm)


# Looping through EPN sample file with both RNA and DNA BSID's   
with open("analyses/molecular-subtyping-EPN/results/EPN_molecular_subtype.tsv", "r") as notebook:
    for line in notebook:
        line = line.split()
        # Adding fusion, CNA, breakpoints density and expression data to header line
        if line[0].startswith("Kids"):
            line.extend(("1q_loss", "9p_loss", "9q_loss", "6p_loss", "6q_loss", "11q_loss", "11q_gain", "NFKB_pathway_GSEAscore", "C11orf95-RELA_fusion", "LTBP3--RELA_fusion", "PTEN--TAS2R1_fusion", "C11orf95--YAP1_fusion", "YAP1--MAMLD1_fusion", "YAP1--FAM118B_fusion", "C11orf95--MAML2_fusion","breaks_density-chromosomal_instability", "RELA_expr_Z-scores", "L1CAM_expr_Zscore","ARL4D_expr_Zscore",  "CLDN1_expr_zscore", "CXorf67_expr_zscore", "TKTL1_expr_zscore", "GPBP1_expr_zscore", "IFT46_expr_zscore"))
            #print(line)
            # Printing header line
            for i in line:
                outfile.write(i+"\t")
            outfile.write("\n")    
        else:
            # Looking for NFKB pathway and fusions based on RNA BSID's 
            nfkb_gsva_score = (gsva_NFKB.get_value(line[3], "gsea_score"))
            C11orf95_RELA = fusion.get_value(line[3], "C11orf95--RELA")
            LTBP3_RELA = fusion.get_value(line[3], "LTBP3--RELA")
            PTEN_TAS2R1 = fusion.get_value(line[3], "PTEN--TAS2R1")
            C11orf95_YAP1 = fusion.get_value(line[3], "C11orf95--YAP1")
            YAP1_MAMLD1 = fusion.get_value(line[3], "YAP1--MAMLD1")
            YAP1_FAM118B = fusion.get_value(line[3], "YAP1--FAM118B")
            C11orf95_MAML2 = fusion.get_value(line[3], "C11orf95--MAML2")
            # Looking for expression data and creating Z-scores  
            RELA_zscore = (fpkm_df.get_value("RELA", line[3])-RELA_mean)/RELA_stdev
            L1CAM_zscore = (fpkm_df.get_value("L1CAM", line[3])-L1CAM_mean)/L1CAM_stdev
            ARL4D_zscore = (fpkm_df.get_value("ARL4D", line[3])-ARL4D_mean)/ARL4D_stdev
            CLDN1_zscore = (fpkm_df.get_value("CLDN1", line[3])-CLDN1_mean)/CLDN1_stdev
            CXorf67_zscore = (fpkm_df.get_value("CXorf67", line[3])-CXorf67_mean)/CXorf67_stdev
            TKTL1_zscore = (fpkm_df.get_value("TKTL1", line[3])-TKTL1_mean)/TKTL1_stdev
            GPBP1_zscore = (fpkm_df.get_value("GPBP1", line[3])-GPBP1_mean)/GPBP1_stdev
            IFT46_zscore = (fpkm_df.get_value("IFT46", line[3])-IFT46_mean)/IFT46_stdev
            # Looking for chromosomal CNA based on DNA BSID's. Some of the corresponding DNA ID's show "NA", so splitting with if statement
            if line[2]!="NA":
                oneqvalue = (CNA.get_value("1q",line[2]))
                nineploss = (CNA.get_value("9p",line[2]))
                nineqloss = (CNA.get_value("9q",line[2]))
                sixploss = (CNA.get_value("6p",line[2]))
                sixqloss = (CNA.get_value("6q",line[2]))
                elevenq = (CNA.get_value("11q",line[2]))
                breaks_density =  breakpoint_density.get_value(line[2], "breaks_density")
                if int(oneqvalue) < 0:
                    line.append("1")
                else:
                    line.append("0")
                if int(nineploss) < 0:
                    line.append("1")
                else:
                    line.append("0")    
                if int(nineqloss) < 0:
                    line.append("1")
                else:
                    line.append("0")    
                if int(sixploss) < 0 :
                    line.append("1")
                else:
                    line.append("0")    
                if int(sixqloss) < 0:
                    line.append("1")
                else:
                    line.append("0")
                if int(elevenq) < 0:
                    line.append("1")
                else:
                    line.append("0")
                if int(elevenq) > 0:
                    line.append("1")
                else:
                    line.append("0")
                # Adding all the data collected so far for every sample with RNA and DNA BSID to "line"     
                line.extend((nfkb_gsva_score, C11orf95_RELA, LTBP3_RELA, PTEN_TAS2R1, C11orf95_YAP1, YAP1_MAMLD1, YAP1_FAM118B, C11orf95_MAML2, breaks_density, RELA_zscore, L1CAM_zscore,ARL4D_zscore, CLDN1_zscore, CXorf67_zscore, TKTL1_zscore, GPBP1_zscore, IFT46_zscore))
                count+=1
                # Writing line to output file 
                for i in line:
                    outfile.write(str(i)+"\t")
                outfile.write("\n")    
            elif line[2]=="NA":
                # Adding all the data collected so far for RNA BSID, All corresponding DNA ID's here are "NA"
                line.extend(("0", "0", "0", "0", "0", "0", "0", nfkb_gsva_score, C11orf95_RELA, LTBP3_RELA, PTEN_TAS2R1, C11orf95_YAP1, YAP1_MAMLD1, YAP1_FAM118B, C11orf95_MAML2, "NA", RELA_zscore, L1CAM_zscore,ARL4D_zscore, CLDN1_zscore, CXorf67_zscore, TKTL1_zscore, GPBP1_zscore, IFT46_zscore))
                count+=1
                #print(line)
                # Writing to output file 
                for i in line:
                    outfile.write(str(i)+"\t")
                outfile.write("\n")
            
print("Done!!")




