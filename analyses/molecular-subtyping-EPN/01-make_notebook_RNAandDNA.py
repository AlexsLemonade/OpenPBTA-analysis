#!/usr/bin/env python
#Author - Teja Koganti (D3B)


import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--histologies', required = True,
                    help = 'path to the histology file')
parser.add_argument('-o', '--outnotebook', required = True,
                    help = "output notebook")
args = parser.parse_args()

pbta_histologies = pd.read_csv(args.histologies, sep="\t")
outnotebook = open(args.outnotebook, "w")

# Based on primary site, supra/infra category  assigned to each sample
# These two primary sites could not be categorized  and will be assigned "None" under disease group
#		"Other locations NOS" and "Ventricles"
def group_disease(primary_site):
    infra = ["posterior fossa",
             "optic",
             "spinal",
             "tectum",
             "spine"]
    supra = ["frontal lobe",
             "parietal lobe",
             "occipital lobe",
             "temporal lobe"]
    primary = primary_site.lower() # this will prevent possible errors from case mismatches
    for site in infra:
        if site in primary:
            return "infratentorial"
    for site in supra:
        if site in primary:
            return "supratentorial"
    # Note we only get to the below return if the primary site was not in either defined group.
    return "undetermined"


# Filtering for ependymoma samples 
EP = pbta_histologies[pbta_histologies["integrated_diagnosis"]=="Ependymoma"]
EP_rnaseq_samples = EP[EP["experimental_strategy"] == "RNA-Seq"][["Kids_First_Biospecimen_ID", "primary_site", 
	"Kids_First_Participant_ID", "sample_id", "experimental_strategy"]]
EP_rnaseq_samples["disease_group"] = [group_disease(primary) for primary in EP_rnaseq_samples["primary_site"]]

# List with only RNA samples
EP_rnasamplenames_PTIDs = list(EP_rnaseq_samples["Kids_First_Participant_ID"]) 

# Filtering for DNA samples 
all_WGS = EP[EP["experimental_strategy"]=="WGS"]
WGSPT = all_WGS[all_WGS["Kids_First_Participant_ID"].isin(EP_rnasamplenames_PTIDs)]
WGS_dnaseqsamples = WGSPT[["Kids_First_Biospecimen_ID", "Kids_First_Participant_ID", "sample_id"]]


# Renaming the column name so they don't conflict in merge step 
EP_rnaseq_samples = EP_rnaseq_samples.rename(columns={"Kids_First_Biospecimen_ID":"Kids_First_Biospecimen_ID_RNA"})
WGS_dnaseqsamples = WGS_dnaseqsamples.rename(columns={"Kids_First_Biospecimen_ID":"Kids_First_Biospecimen_ID_DNA"})


# sample_id is common between both  datafarmes and also unique between RNA and DNA. 
# Some DNA BSID's are missing for the corresponding RNA samples 
EP_rnaseq_WGS = EP_rnaseq_samples.merge(WGS_dnaseqsamples, 
                                        on = ["sample_id", "Kids_First_Participant_ID"], 
                                        how = "outer")
EP_rnaseq_WGS.fillna('NA', inplace=True)

EP_rnaseq_WGS[["Kids_First_Participant_ID", 
		"sample_id", 
		"Kids_First_Biospecimen_ID_DNA", 
		"Kids_First_Biospecimen_ID_RNA", 
		"disease_group"
		]].to_csv(outnotebook, sep="\t", index=False)
outnotebook.close()
    
