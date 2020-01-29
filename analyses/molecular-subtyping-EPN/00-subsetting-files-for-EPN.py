#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np 



pip install rpy2



import  pyreadr
import rpy2.robjects as robjects


#Reading in pbta  histologies file to subset just  the ependymoma samples
pbta_histologies = pd.read_csv("data/pbta-histologies.tsv", sep="\t")


EP = pbta_histologies[pbta_histologies["disease_type_new"]=="Ependymoma"]
EP_rnaseq_samples = EP[EP["experimental_strategy"] == "RNA-Seq"][["Kids_First_Biospecimen_ID", "primary_site"]]
EP_rnaseq_samples["disease_group"] = ["infratentorial" if "Posterior Fossa" in primary else "infratentorial" if "Optic" in primary else "supratentorial" if "Frontal Lobe" in primary else "supratentorial" if "Parietal Lobe" in primary else "infratentorial" if "Spinal Cord" in primary else "supratentorial" if "Occipital Lobe" in primary else "infratentorial" if "Tectum" in primary else "infratentorial" if "Spine" in primary else "supratentorial" if "Temporal Lobe" in primary else "infratentorial" if "Spinal" in primary else  "None" for primary in EP_rnaseq_samples["primary_site"]]

# This list only has BSID;s for ependymoma samples
EP_rnasamplenames = list(EP_rnaseq_samples["Kids_First_Biospecimen_ID"])


import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()

# Reading in pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds file to subset (All ependymoma samples are  stranded, so ignoring polyA gene expression file in subsetting )
readRDS = robjects.r['readRDS']
collapsed_rsem_stranded = readRDS("data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds")

# Subsetting columns with column names/BSIDs that are  in the  list  of ependymoma samples
collapsed_rsem_stranded_only_with_ependymomasamples = collapsed_rsem_stranded.reindex(columns = EP_rnasamplenames)

saveRDS = robjects.r['saveRDS']
saveRDS(collapsed_rsem_stranded_only_with_ependymomasamples, "analyses/molecular-subtyping-EPN/epn-subset/epn_pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds")



# In[156]:





# In[158]:





