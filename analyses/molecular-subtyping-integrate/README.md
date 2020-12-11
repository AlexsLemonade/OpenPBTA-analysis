## Integrate molecular subtyping output from pathology feedback

**Author of code and documentation:** [@kgaonkar6](https://github.com/kgaonkar6)

In this repo we add molecular subtype for various broad_histologies but in some cases, results are updated in `molecular-subtyping-pathology`. In this analysis we add molecular_subtype,	integrated_diagnosis, short_histology,	broad_histology and	Notes from  `compiled_molecular_subtypes_with_clinical_pathology_feedback.tsv`.

### Usage
```sh
bash run-subtyping-integrate.sh
```

### Module contents

`01-integrate-subtyping.Rmd` integrates results from compiled results in `compiled_molecular_subtypes_with_clinical_pathology_feedback.tsv` to `pbta-histologies-base.tsv`

