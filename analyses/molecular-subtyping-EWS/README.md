## Molecular Subtyping of Ewing's Sarcoma  

**Module authors:** Krutika Gaonkar ([@kgaonkar6](https://github.com/kgaonkar6)), Jo Lynne Rokita ([@jharenza](https://github.com/jharenza)) and Jaclyn Taroni ([@jaclyn-taroni](https://github.com/jaclyn-taroni))

_EWSR1_ fusions are hallmark alterations of Ewing's Sarcoma but we see these fusions in samples with other broad_histologies in our dataset (Diffuse astrocytic and oligodendroglial tumor (1), Metastatic secondary tumors (2), Tumors of sellar region (1)). 
This analysis will subtype these samples as "EWS" (Ewing's Sarcoma) tumors.

### Usage

This module can be run with the following:

```
bash run_subtyping.sh
```

### Module contents

`01-run-subtyping-ewings.Rmd` identifies sample_id for RNA-Seq samples with hallmark _EWSR1_ fusions and matches sample_id with WGS sample and adds molecular_subtype as "EWS" (Ewing Sarcoma) 
The tabular output (`results/EWS_results.tsv`) only contains tumors where `molecular_subtype == "EWS"`.
