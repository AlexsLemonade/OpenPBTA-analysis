#!/usr/bin/env python


"""
01-snv-frequencies.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Functions to create copy number variation (CNV) cancer type and study gene-level frequencies for OPenPedCan analyses modules 
"""


__author__ = ('Eric Wafula (wafulae@chop.edu)')
__version__ = '1.0'
__date__ = '12 July 2021'


import os 
import sys
import csv
import json
import argparse
import numpy as np
import pandas as pd
from functools import reduce
from collections import OrderedDict


def read_parameters():
     p = argparse.ArgumentParser(description=("The 01-snv-frequencies.py scripts creates copy number variation (CNV) cancer type and study gene-level alterations frequencies table for the OPenPedCan analyses modules."), formatter_class=argparse.RawTextHelpFormatter)
     p.add_argument('HISTOLOGY_FILE', type=str, default=None, help="OPenPedCan histology file (histologies.tsv)\n\n")
     p.add_argument('CNV_FILE', type=str, default=None, help="OPenPedCan CNV consensus file (consensus_seg_annotated_cn_autosomes.tsv or consensus_seg_annotated_cn_x_and_y.tsv)\n\n")
     p.add_argument('PRIMARY_TUMORS', type=str, default=None, help="OPenPedCan independent primary tumor samples file (independent-specimens.wgs.primary.tsv)\n\n")
     p.add_argument('RELAPSE_TUMORS', type=str, default=None, help="OPenPedCan independent relapse tumor samples file (independent-specimens.wgs.relapse.tsv)\n\n")
     p.add_argument('ENSG_NAME', type=str, default=None, help="OPenPedCan Ensembl to gene full name mapping file (ensg-gene-full-name-refseq-protein.tsv)\n\n")
     p.add_argument('ONCOKB', type=str, default=None, help="OPenPedCan Hugo gene symbols to OncoKB categories mapping file (OncoKB_Oncogene-TSG_genes.tsv)\n\n")
     p.add_argument('EFO_MONDO', type=str, default=None, help="OPenPedCan disease to EFO/MONDO mapping file (efo-mondo-map.tsv)\n\n")
     p.add_argument('ENSG_RMTL', type=str, default=None, help="OPenPedCan Ensembl to RMTL mapping file (ensg-hugo-rmtl-v*-mapping.tsv)\n\n")
     p.add_argument('-v', '--version', action='version', version="classify_reads.py version {} ({})".format(__version__, __date__), help="Print the current 01-cnv-frequencies.py version and exit\n\n")
     return p.parse_args()


def merge_histology_and_cnv_data(histology_file, cnv_consensus_file):
     # load histology file
     histology_df = pd.read_csv(histology_file, sep="\t", na_filter=False, dtype=str)
     
     # subset histology dataframe for relevant columns
     histology_df = histology_df[["Kids_First_Biospecimen_ID","Kids_First_Participant_ID", "cohort", "cancer_group", "sample_type"]]
     
     # load CNV consensus file
     cnv_df = pd.read_csv(cnv_consensus_file, sep="\t")
     
     # merge subset of histology dataframe to CNV dataframe keeping only sample present in the CNV table (left outer join)
     merged_df = pd.merge(cnv_df, histology_df, how="left", left_on="biospecimen_id", right_on="Kids_First_Biospecimen_ID")
     
     # check if non tumor sample are present in the merged dataframe
     if not np.array_equal(np.array(["Tumor"]), merged_df.sample_type.unique()):
          raise Exception("Merged hsitology-CNV dataframe contains non tumor samples")
     
     # check and drop unknown cancer types (cancer_group == NA)
     row_indices = merged_df[(merged_df["cancer_group"] == "NA")].index
     merged_df.drop(row_indices, inplace=True)
     
     # select and reorder relevant columns
     all_tumors_df = merged_df[["gene_symbol", "ensembl", "status", "copy_number", "ploidy", "cohort", "cancer_group", "Kids_First_Biospecimen_ID", "Kids_First_Participant_ID"]].reset_index()
     return(all_tumors_df)


def get_cancer_groups_and_cohorts(all_tumors_df):
     # group samples by cancer_group and cohort
     cancer_group_cohort_df = all_tumors_df.groupby(["cancer_group", "cohort"])["Kids_First_Biospecimen_ID"].nunique().reset_index()
     cancer_group_cohort_df.columns = ["cancer_group", "cohort", "num_samples"]
     
     # group samples by cancer_group only
     def func(x):
          d = {}
          cohort_list = x["cohort"].unique()
          if len(cohort_list) > 1:
               d["cohort"] = "all_cohorts"
          else:
               d["cohort"] = cohort_list[0]
          d["num_samples"] = x["Kids_First_Biospecimen_ID"].nunique()
          return pd.Series(d, index=["cohort", "num_samples"])
     cancer_group_df = all_tumors_df.groupby(["cancer_group"]).apply(func).reset_index()
     
     # concat cancer_group_cohort_df and cancer_group_df and rremove duplicates
     cancer_group_cohort_df = pd.concat([cancer_group_df, cancer_group_cohort_df], sort=False, ignore_index=True)
     cancer_group_cohort_df.drop_duplicates(inplace=True)
     # these subset dataframe is for testing and need to comment out
     #cancer_group_cohort_df = cancer_group_cohort_df[cancer_group_cohort_df.cancer_group.isin(["CNS Embryonal tumor", "Atypical Teratoid Rhabdoid Tumor"])]
     return(cancer_group_cohort_df)


def compute_variant_frequencies(all_tumors_df, primary_tumors_file, relapase_tumors_file, cancer_group_cohort_df):
     # get independent primary tumor samples and subset from all tumor sample dataframe
     tumor_dfs = {"all_tumors": all_tumors_df}
     primary_tumors_samples_list = list(pd.read_csv(primary_tumors_file, sep="\t")["Kids_First_Biospecimen_ID"].unique())
     primary_tumors_df = all_tumors_df[all_tumors_df.Kids_First_Biospecimen_ID.isin(primary_tumors_samples_list)].reset_index()
     tumor_dfs["primary_tumors"] = primary_tumors_df
     
     # get independent relapse tumor sample and subset from all tumor samples dataframe
     relapse_tumors_samples_list = list(pd.read_csv(relapase_tumors_file, sep="\t")["Kids_First_Biospecimen_ID"].unique())
     relapse_tumors_df = all_tumors_df[all_tumors_df.Kids_First_Biospecimen_ID.isin(relapse_tumors_samples_list)].reset_index()
     tumor_dfs["relapse_tumors"] = relapse_tumors_df
     
     # compute variant frequencies for each cancer group per cohort and  cancer group in cohorts
     # for the overal dataset (all tumor samples)  and independent primary/replase tumor samples
     def func(x):
          d = {}
          d["Gene_Symbol"] = ",".join(x["gene_symbol"].unique())
          d["total_sample_alterations"] = x["Kids_First_Biospecimen_ID"].nunique()
          d["total_patient_alterations"] = x["Kids_First_Participant_ID"].nunique()
          return(pd.Series(d, index=["Gene_Symbol", "total_sample_alterations", "total_patient_alterations"]))
     all_tumors_frequecy_dfs = []
     primary_tumors_frequecy_dfs = []
     relapse_tumors_frequecy_dfs = []
     for index, row in cancer_group_cohort_df.iterrows():
          if row["num_samples"] > 5:
               for df_name, tumor_df in tumor_dfs.items():
                    if row["cohort"] == "all_cohorts":
                         df = tumor_df[(tumor_df["cancer_group"] == row["cancer_group"])]
                    else:
                         df = tumor_df[(tumor_df["cancer_group"] == row["cancer_group"]) & (tumor_df["cohort"] == row["cohort"])]
                    if df.empty:
                         continue
                    num_samples = df["Kids_First_Biospecimen_ID"].nunique()
                    num_patients = df["Kids_First_Participant_ID"].nunique()
                    df = df.groupby(["ensembl", "status"]).apply(func)
                    df = df.rename_axis(["Gene_Ensembl_ID", "Variant_Type"]).reset_index()
                    df["num_patients"] = num_patients
                    for i, j in df.iterrows():
                         if df_name == "primary_tumors" or df_name == "relapse_tumors":
                              df.at[i, "Total_Alterations/Patients_in_Dataset"] = "{}/{}".format(j["total_sample_alterations"], num_samples)
                              df.at[i, "Frequency_in_Overall_Dataset"] = "{:.2f}%".format((j["total_sample_alterations"]/num_samples)*100)
                         if df_name == "all_tumors":
                              df.at[i, "Total_Alterations/Patients_in_Dataset"] = "{}/{}".format(j["total_patient_alterations"], num_patients)
                              df.at[i, "Frequency_in_Overall_Dataset"] = "{:.2f}%".format((j["total_patient_alterations"]/num_patients)*100)
                         df.at[i, "Dataset"] = row["cohort"]
                         df.at[i, "Disease"] = row["cancer_group"]
                    df = df [["Gene_Symbol", "Gene_Ensembl_ID", "Variant_Type", "Dataset", "Disease", "Total_Alterations/Patients_in_Dataset", "Frequency_in_Overall_Dataset"]]
                    if df_name == "primary_tumors":
                         primary_tumors_frequecy_dfs.append(df)
                    elif df_name == "relapse_tumors":
                         relapse_tumors_frequecy_dfs.append(df)
                    else:
                         all_tumors_frequecy_dfs.append(df)
                         
     # merge overal dataset (all tumor samples) and independent primary/replase tumor samples
     # frequencies for cancer groups per cohorts into a single dataframe
     merging_list = []
     #frequencies in overall dataset
     all_tumors_frequecy_df = pd.concat(all_tumors_frequecy_dfs, sort=False, ignore_index=True)
     merging_list.append(all_tumors_frequecy_df)
     # frequencies in independent primary tumors
     primary_tumors_frequecy_df = pd.concat(primary_tumors_frequecy_dfs, sort=False, ignore_index=True)
     primary_tumors_frequecy_df.rename(columns={"Total_Alterations/Patients_in_Dataset": "Total_Primary_Tumors_Altered/Primary_Tumors_in_Dataset", "Frequency_in_Overall_Dataset": "Frequency_in_Primary_Tumors"}, inplace=True)
     merging_list.append(primary_tumors_frequecy_df)
     # frequencies in independent relapse tumors
     relapse_tumors_frequecy_df = pd.concat(relapse_tumors_frequecy_dfs, sort=False, ignore_index=True)
     relapse_tumors_frequecy_df.rename(columns={"Total_Alterations/Patients_in_Dataset": "Total_Relapse_Tumors_Altered/Relapse_Tumors_in_Dataset", "Frequency_in_Overall_Dataset": "Frequency_in_Relapse_Tumors"}, inplace=True)
     merging_list.append(relapse_tumors_frequecy_df)
     cnv_frequency_df = reduce(lambda x, y: pd.merge(x, y, how="outer", on=["Gene_Symbol", "Gene_Ensembl_ID", "Variant_Type", "Dataset", "Disease"]), merging_list).fillna("")
     cnv_frequency_df = cnv_frequency_df.replace({"Total_Primary_Tumors_Altered/Primary_Tumors_in_Dataset": "", "Total_Relapse_Tumors_Altered/Relapse_Tumors_in_Dataset": ""}, "0/0")
     return(cnv_frequency_df)


def get_annotations(cnv_frequency_df, ensg_gene_name_mapping_file, oncokb_mapping_file, efo_mondo_mapping_file, rmtl_mapping_file):
     # insert annotation columns in the CNV frequency dataframe
     cnv_frequency_df["RMTL"] = cnv_frequency_df.insert(1, "RMTL", "")
     cnv_frequency_df["Gene_Full_Name"] = cnv_frequency_df.insert(3, "Gene_Full_Name", "")
     cnv_frequency_df["Variant_Category"] = cnv_frequency_df.insert(5, "Variant_Category", "")
     cnv_frequency_df["OncoKB_Category"] = cnv_frequency_df.insert(6, "OncoKB_Category", "")
     cnv_frequency_df["EFO_ID"] = cnv_frequency_df.insert(7, "EFO_ID", "")
     cnv_frequency_df["MONDO_ID"] = cnv_frequency_df.insert(8, "MONDO_ID", "")
     cnv_frequency_df.fillna("", inplace=True)

     # populate CNV frequency dataframe with full gene names
     row_index = -1
     ensembl_gene_name_df = pd.read_csv(ensg_gene_name_mapping_file, sep="\t")
     for ensembl_id in ensembl_gene_name_df["Gene_Ensembl_ID"]:
          row_index += 1
          if ensembl_id in list(cnv_frequency_df.Gene_Ensembl_ID):
               cnv_frequency_df.loc[cnv_frequency_df.Gene_Ensembl_ID == ensembl_id, "Gene_Full_Name"] = ensembl_gene_name_df.at[row_index, "Gene_full_name"]

     # populate CNV frequency dataframe with OncoKB categories
     row_index = -1
     oncokb_df = pd.read_csv(oncokb_mapping_file, sep="\t")
     for gene_symbol in oncokb_df["hugo_symbol"]:
          row_index += 1
          if gene_symbol in list(cnv_frequency_df.Gene_Symbol):
               cnv_frequency_df.loc[cnv_frequency_df.Gene_Symbol == gene_symbol, "OncoKB_Category"] = oncokb_df.at[row_index, "oncokb_category"]

     # populate CNV frequency dataframe with EFO and MONDO disease accessions
     row_index = -1
     efo_mondo_df = pd.read_csv(efo_mondo_mapping_file, sep="\t")
     for cancer_group in efo_mondo_df["cancer_group"]:
          row_index += 1
          if cancer_group in list(cnv_frequency_df.Disease):
               cnv_frequency_df.loc[cnv_frequency_df.Disease == cancer_group , "EFO_ID"] = efo_mondo_df.at[row_index, "efo_code"]
               cnv_frequency_df.loc[cnv_frequency_df.Disease == cancer_group , "MONDO_ID"] = efo_mondo_df.at[row_index, "mondo_code"]

     # populate CNV frequency dataframe with EFO and MONDO disease accessions
     row_index = -1
     ensg_rmtl_df = pd.read_csv(rmtl_mapping_file, sep="\t")
     ensg_rmtl_df = ensg_rmtl_df[ensg_rmtl_df.rmtl.isin(["Relevant Molecular Target"])].reset_index()
     for ensembl_gene in ensg_rmtl_df["ensg_id"]:
          row_index += 1
          if ensembl_gene in list(cnv_frequency_df.Gene_Ensembl_ID):
               cnv_frequency_df.loc[cnv_frequency_df.Gene_Ensembl_ID == ensembl_gene , "RMTL"] = ensg_rmtl_df.at[row_index, "rmtl"]
     cnv_frequency_df.fillna("", inplace=True)
     return(cnv_frequency_df)


def main():
     # get input parameters
     args = read_parameters()

     # call functions to compute CNV gene-level and add functional annotations
     all_tumors_df = merge_histology_and_cnv_data(args.HISTOLOGY_FILE, args.CNV_FILE)
     cancer_group_cohort_df = get_cancer_groups_and_cohorts(all_tumors_df)
     cnv_frequency_df = compute_variant_frequencies(all_tumors_df, args.PRIMARY_TUMORS, args.RELAPSE_TUMORS, cancer_group_cohort_df)
     cnv_frequency_df = get_annotations(cnv_frequency_df, args.ENSG_NAME, args.ONCOKB, args.EFO_MONDO, args.ENSG_RMTL)

     # write annotated CNV frequencies to the results directory
     results_dir = "{}/results".format(os.path.dirname(__file__))
     if not os.path.exists(results_dir):
          os.mkdir(results_dir)
     # write results to TSV file
     results_tsv = "{}/{}_freq.txt".format(results_dir, os.path.splitext(os.path.splitext(os.path.basename(args.CNV_FILE))[0])[0])
     cnv_frequency_df.to_csv(results_tsv, sep="\t", index=False, encoding="utf-8")
     # write results to TSV file
     results_jsonl = "{}/{}_freq.jsonl".format(results_dir, os.path.splitext(os.path.splitext(os.path.basename(args.CNV_FILE))[0])[0])
     tsv_file = open(results_tsv, 'r')
     jsonl_file = open(results_jsonl, 'w')
     reader = csv.DictReader(tsv_file, delimiter="\t")
     headers = reader.fieldnames
     for row in reader:
          row_dict = OrderedDict()
          for header in headers:
               row_dict[header] = row[header]
          json.dump(row_dict, jsonl_file)
          jsonl_file.write("\n")
     tsv_file.close()
     jsonl_file.close()
     sys.exit(0)


if __name__ == '__main__':
     main()	
