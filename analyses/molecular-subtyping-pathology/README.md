## Updating molecular subtyping output with pathology feedback

**Author of code and documentation:** [@jaclyn-taroni](https://github.com/jaclyn-taroni)

As part of this project, we have undertaken several analyses that use molecular data to subtype biospecimens. 
In some cases, our analyses have resulted in an update of the `integrated_diagnosis` field included in the clinical file (`pbta-histologies.tsv` [[doc](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/data-formats.md#data-caveats)]).

The objective of this module is three-fold:

1. Aggregate the molecular subtyping calls from the following modules that produced results:
   * [`molecular-subtyping-EPN`](https://github.com/jaclyn-taroni/OpenPBTA-analysis/tree/645-pathology-feedback/analyses/molecular-subtyping-EPN)
   * [`molecular-subtyping-EWS`](https://github.com/jaclyn-taroni/OpenPBTA-analysis/tree/645-pathology-feedback/analyses/molecular-subtyping-EWS)
   * [`molecular-subtyping-HGG`](https://github.com/jaclyn-taroni/OpenPBTA-analysis/tree/645-pathology-feedback/analyses/molecular-subtyping-HGG)
   * [`molecular-subtyping-LGAT`](https://github.com/jaclyn-taroni/OpenPBTA-analysis/tree/645-pathology-feedback/analyses/molecular-subtyping-LGAT)
   * [`molecular-subtyping-embryonal`](https://github.com/jaclyn-taroni/OpenPBTA-analysis/tree/645-pathology-feedback/analyses/molecular-subtyping-embryonal)
   * [`molecular-subtyping-CRANIO`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/molecular-subtyping-CRANIO)
   * [`molecular-subtyping-MB`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/molecular-subtyping-MB)
   * [`molecular-subtyping-neurocytoma`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/molecular-subtyping-neurocytoma)
   * [`molecular-subtyping-embryonal`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/molecular-subtyping-embryonal)
   * [`molecular-subtyping-chordoma`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/molecular-subtyping-chordoma)

2. Incorporate clinical reviewed subtypes for PNOC003 samples:
In the original [issue](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/751) added by @jharenza we have the clinically reviewed subtypes for PNOC003 samples. We check if any subtype is different between the results from `molecular-subtyping-HGG` and this file and update to the clinically reviewed subtype. Subtypes for 3 WXS samples and 3 RNA-Seq from PT_NK8A49X5, PT_QA9WJ679 and PT_WGVEF96B were updated.  

3. Incorporate feedback from CHOP pathologists Maria Rita Santi and Angela Viaene. 
Specifically, there are instances where the final `integrated_diagnosis` calls from pathology will deviate from the logic included in molecular subtyping modules based on additional information outside the scope of the repository (e.g., pathology reports, slides, etc.). 
The goal is to make sure that the _final calls_ are recorded in an aggregated table (see point 1 above) and documented in this repository.

For more background, see [#609](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/609).

Chordoma samples needed pathology review to add molecular subtypes at this step as well, see [#608](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/608).  


### Usage

To run all steps in this module, run the following command:

```sh
bash run-subtyping-aggregation.sh
```

### Module contents

`01-compile-subtyping-results.Rmd` aggregates results from the modules listed above into a single table (`results/compiled_molecular_subtypes.tsv`).

`02-incorporate-clinical-feedback.Rmd` incorporate clincally reviewed subtypes for PNOC003 samples and update to the clinically reviewed subtype if they are different from the subtype from `molecular-subtyping-HGG`

`03-incorporate-pathology-feedback.Rmd` incorporates pathology feedback for specific samples when the labels for those samples either need to be updated as a result of molecular subtyping OR molecular abberations data could not idenitify subtypes OR pathology review deviates from the logic in upstream molecular subtyping modules. The output is an updated version of the table from `01-compile-subtyping-results.Rmd` (`results/compiled_molecular_subtypes_with_clinical_pathology_feedback.tsv`).
