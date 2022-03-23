## Integrate molecular subtyping output from pathology feedback

**Author of code and documentation:** [@kgaonkar6](https://github.com/kgaonkar6)

In this repo, we add molecular subtype for molecular_subtype from all subtyping modules and integrated_diagnosis, short_histology, broad_histology, and Notes from `compiled_molecular_subtypes_with_clinical_pathology_feedback.tsv`. A column `cancer_group` is added to provide broader terms derived from `harmonized_diagnosis` which can be used to generate figures.

### Usage
```sh
bash run-subtyping-integrate.sh
```

### Module contents

`01-integrate-subtyping.Rmd` integrates results from compiled results in `compiled_molecular_subtypes_with_clinical_pathology_feedback.tsv` to `pbta-histologies-base.tsv`
To assign `cancer_group`, we follow a two-step procedure:
1) remove molecular subtype information from `harmonized_diagnosis`
2) make additional modifications based on the subtype-stripped diagnosis (largely formatting changes)

Additional modifications after subtype removal are as followed:
| cancer group pre-modification (with or without subtype removal) | final cancer group used | 
|-----------|----------------|
| Adamantinomatous craniopharyngioma | Craniopharyngioma |
| Anaplastic (malignant) meningioma | Meningioma |
| Atypical meningioma | Meningioma |
| Atypical Teratoid Rhabdoid Tumor (ATRT) | Atypical Teratoid Rhabdoid Tumor |
| Benign tumor | NA |
| Brainstem glioma- Diffuse intrinsic pontine glioma | Diffuse intrinsic pontine glioma |
| Clear cell meningioma | Meningioma |
| CNS Embryonal tumor | Embryonal tumor with multilayer rosettes |
| Dysembryoplastic neuroepithelial tumor (DNET) | Dysembryoplastic neuroepithelial tumor |
| Dysembryoplastic neuroepithelial tumor (DNET);Dysplasia/Gliosis;Ganglioglioma | Dysembryoplastic neuroepithelial tumor |
| Dysembryoplastic neuroepithelial tumor (DNET);Ganglioglioma | Dysembryoplastic neuroepithelial tumor |
| Dysplasia/Gliosis | NA |
| Dysplasia/Gliosis;Glial-neuronal tumor NOS | Dysplasia Gliosis-Glial-neuronal tumor NOS |
| Germinoma;Teratoma | Germinoma-Teratoma |
| High-grade glioma/astrocytoma | High-grade glioma astrocytoma |
| High-grade glioma/astrocytoma (WHO grade III/IV) | High-grade glioma astrocytoma |
| Low-grade glioma/astrocytoma | Low-grade glioma astrocytoma |
| Low-grade glioma/astrocytoma (WHO grade I/II) | Low-grade glioma astrocytoma |
| Malignant peripheral nerve sheath tumor (MPNST) | Malignant peripheral nerve sheath tumor |
| Meningothelial meningioma | Meningioma |
| Metastatic secondary tumors;Neuroblastoma | Metastatic secondary tumors-Neuroblastoma |
| Neurofibroma/Plexiform | Neurofibroma Plexiform |
| Neurofibroma/Plexiform;Other | Neurofibroma Plexiform |
| Non-germinomatous germ cell tumor;Teratoma | Teratoma | 

Additionally, additional tsv files indicating how cancer groups were matched were generated:
1) `cancer_group_post_subtype_removal_match.tsv`: this file showed the mapping between harmonized diagnosis WITHOUT subtypes (after step 1) and final cancer group.
2) `cancer_group_harmonized_diag_match.tsv`: this file showed the mapping between harmonized diagnosis (before step 1) and final cancer group.

