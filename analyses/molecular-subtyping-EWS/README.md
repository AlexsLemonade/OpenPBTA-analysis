# Molecular Classifying Ewings Sarcoma  

**Module authors:** Krutika Gaonkar ([@kgaonkar6](https://github.com/kgaonkar6)) , Jo Lynne Rokita ([@jharenza](https://github.com/jharenza)) and Jaclyn Taroni ([@jaclyn-taroni](https://github.com/jaclyn-taroni))

_EWSR1_ fusions are hallmark alterations of Ewing Sarcoma but we see these fusions in samples with other broad_histologies in our dataset (Diffuse astrocytic and oligodendroglial tumor (1), Metastatic secondary tumors (2), Tumors of sellar region (1)). This analysis will reclassify these samples as Ewings Sarcoma tumors.


`01-reclassify_as_ewings.Rmd` Identifies `integrated_diagnosis`, `short_histology` from clinical file for sample IDs with hallmark _EWSR1_ fusions and reclassifies as Ewing Sarcoma (integrated_diagnosis) and EWS (short_histology)
