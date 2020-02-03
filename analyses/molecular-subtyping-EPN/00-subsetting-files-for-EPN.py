#!/usr/bin/env python
# coding: utf-8
# Author - Teja Koganti (D3B)

import pandas as pd
import pyreadr


#Reading in pbta  histologies file to subset just  the ependymoma samples
pbta_histologies = pd.read_csv("data/pbta-histologies.tsv", sep="\t")


EP = pbta_histologies[pbta_histologies["disease_type_new"]=="Ependymoma"]
EP_rnaseq_samples = EP[EP["experimental_strategy"] == "RNA-Seq"][["Kids_First_Biospecimen_ID", "primary_site"]]
EP_rnaseq_samples["disease_group"] = ["infratentorial" if "Posterior Fossa" in primary else "infratentorial" if "Optic" in primary else "supratentorial" if "Frontal Lobe" in primary else "supratentorial" if "Parietal Lobe" in primary else "infratentorial" if "Spinal Cord" in primary else "supratentorial" if "Occipital Lobe" in primary else "infratentorial" if "Tectum" in primary else "infratentorial" if "Spine" in primary else "supratentorial" if "Temporal Lobe" in primary else "infratentorial" if "Spinal" in primary else  "None" for primary in EP_rnaseq_samples["primary_site"]]

# This list only has BSID;s for ependymoma samples
EP_rnasamplenames = list(EP_rnaseq_samples["Kids_First_Biospecimen_ID"])

# Reading in pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds file to subset (All ependymoma samples are  stranded, so ignoring polyA gene expression file in subsetting )
collapsed_rsem_stranded = pyreadr.read_r("data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds")[None]

# Subsetting columns with column names/BSIDs that are  in the  list  of ependymoma samples
collapsed_rsem_stranded_only_with_ependymomasamples = collapsed_rsem_stranded.reindex(columns = EP_rnasamplenames)

pyreadr.write_rds("analyses/molecular-subtyping-EPN/epn-subset/epn_pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds", collapsed_rsem_stranded_only_with_ependymomasamples)

