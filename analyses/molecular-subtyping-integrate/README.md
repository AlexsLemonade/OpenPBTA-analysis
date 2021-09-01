## Integrate molecular subtyping output from pathology feedback

**Author of code and documentation:** [@kgaonkar6](https://github.com/kgaonkar6)

In this repo, we add molecular subtype for molecular_subtype from all subtyping modules and integrated_diagnosis, short_histology, broad_histology, and Notes from `compiled_molecular_subtypes_with_clinical_pathology_feedback.tsv`. A column `cancer_group` is added to provide broader terms derived from `harmonized_diagnosis` which can be used to generate figures.

### Usage
```sh
bash run-subtyping-integrate.sh
```

### Module contents

`01-integrate-subtyping.Rmd` integrates results from compiled results in `compiled_molecular_subtypes_with_clinical_pathology_feedback.tsv` to `pbta-histologies-base.tsv`
