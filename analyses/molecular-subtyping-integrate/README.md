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
2) did additional modifications based on the subtype-stripped diagnosis

1. The following `cancer_group` do not have molecular subtypes information and hence, `cancer_group` is first identified as `harmonized_diagnosis` in step 1).
These `cancer_group` (without additional modification on step 2) yet) are:
```
Adamantinomatous craniopharyngioma
Adenoma
Anaplastic (malignant) meningioma
Arteriovenous malformation
Atypical choroid plexus papilloma
Atypical meningioma
Atypical Teratoid Rhabdoid Tumor (ATRT)
Brainstem glioma- Diffuse intrinsic pontine glioma
Cavernoma
Choroid plexus carcinoma
Choroid plexus cyst
Choroid plexus papilloma
Clear cell meningioma
CNS Burkitt’s lymphoma
CNS neuroblastoma
Craniopharyngioma
Diffuse fibrillary astrocytoma
Dysembryoplastic neuroepithelial tumor (DNET)
Dysembryoplastic neuroepithelial tumor (DNET);Dysplasia/Gliosis;Ganglioglioma
Dysembryoplastic neuroepithelial tumor (DNET);Ganglioglioma
Dysplasia/Gliosis
Dysplasia/Gliosis;Glial-neuronal tumor NOS
Ependymoma
Ewing sarcoma
Fibromyxoid lesion
Ganglioglioma
Ganglioneuroblastoma
Ganglioneuroma
Germinoma
Germinoma;Teratoma
Glial-neuronal tumor NOS
Hemangioblastoma
High-grade glioma/astrocytoma (WHO grade III/IV)
Juvenile xanthogranuloma
Langerhans Cell histiocytosis
Low-grade glioma/astrocytoma (WHO grade I/II)
Malignant peripheral nerve sheath tumor (MPNST)
Medulloblastoma
Melanocytic tumor
Meningioma
Meningothelial meningioma
Metastatic secondary tumors
Metastatic secondary tumors;Neuroblastoma
Myofibroblastoma
Myxoid spindle cell tumor
Neuroblastoma
Neurocytoma
Neurofibroma/Plexiform
Neurofibroma/Plexiform;Other
Non-germinomatous germ cell tumor;Teratoma
Oligodendroglioma
Pilocytic astrocytoma
Pineoblastoma
Pleomorphic xanthoastrocytoma
Reactive connective tissue
Rhabdomyosarcoma
Rosai-Dorfman disease
Sarcoma
Schwannoma
Subependymal Giant Cell Astrocytoma
Teratoma
```
2. The following `cancer_group` have molecular subtypes information and hence, `cancer_group` is first identified as `harmonized_diagnosis` without subtype after step 1).
These `cancer_group` (after molecular subtype removal without additional modification from step 2) yet) are:
```
Benign tumor
Chordoma
CNS Embryonal tumor
Desmoplastic infantile astrocytoma and ganglioglioma
Diffuse leptomeningeal glioneuronal tumor
Diffuse midline glioma
Embryonal tumor with multilayer rosettes
Ependymoma
Epilepsy
Ganglioglioma
Glial-neuronal tumor NOS
High-grade glioma/astrocytoma
Low-grade glioma/astrocytoma
Medulloblastoma
Pleomorphic xanthoastrocytoma
Rosette-forming glioneuronal tumor
Subependymal Giant Cell Astrocytoma
```

3. After step 1), we made additional modifications as step 2).
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
