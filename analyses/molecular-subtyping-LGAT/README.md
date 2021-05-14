# Molecular Subtyping of Low Grade Gliomaa

**Module authors:** Krutika Gaonkar ([@kgaonkar6](https://github.com/kgaonkar6)), Jaclyn Taroni ([@jaclyn-taroni](https://github.com/jaclyn-taroni))

In this analysis we subtype LGAT samples according to the presence/absence of molecular alterations decribed below and originally discussed in #790(https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/790):

- **LGG, NF1-somatic**
  - This subtype is characterized by somatic _NF1_ variants.

- **LGG, NF1-germline**
  - We do not have germline data within OpenPBTA, but we do have information about cancer predispositions. This group can be annotated by checking for `NF-1` within `cancer_predispositions` field of the `pbta-histologies.tsv` file. 

- **LGG, KIAA1549-BRAF**
  - contains _KIAA1549-BRAF_ fusion

- **LGG, BRAF V600E**
  - contains BRAF V600E or V599 SNV or non-canonical _BRAF_ alterations such as p.V600ins or p.D594N within the kiinase (PK_Tyr_Ser-Thr) domain.

- **LGG, other MAPK**
  - contains non-canonical _BRAF_ fusion other than _KIAA1549-BRAF_
  - contains _RAF1_ fusion
  - contains _KRAS_, _NRAS_, _HRAS_, _MAP2K1_, _MAP2K2_, _MAP2K1_, _ARAF_, _RAF1_ and BRAF (other than V600E or V599 SNV and other mutations in kinase domain) SNV or indel

- **LGG, RTK**
  - harbors a fusion in _ALK_, _ROS1_, _NTRK1_, _NTRK2_, or _NTRK3_ **or**
  - harbors a _MET_ SNV **or 
  - harbors a _KIT_ SNV **or**
  - harbors a _PDGFRA_ SNV or fusion

- **LGG, FGFR**
  - harbors _FGFR1_ p.N546K, p.K656E, p.N577, or p. K687 hotspot mutations **or** 
  - harbors _FGFR1_ TKD (tyrosine kinase domain tandem duplication) **or**
  - harbors _FGFR1-TACC1_ fusion **or 
  - harbors _FGFR1_ or _FGFR2_ fusions

- **LGG, IDH**
  - harbors an _IDH_ R132 mutation

- **LGG, H3**
  - harbors an _H3F3A_ K28M or G35R/V mutation **or**
  - harbors an _H3F3B_,_HIST1H3B_,_HIST1H3C_,_HIST2H3C_ K28M mutation

- **LGG, MYB/MYBL1**
  - harbors either a _MYB-QKI_ fusion or other _MYB_ or _MYBL1_ fusion

- **LGG, CDKN2A/B** 
  - harbors focal CDKN2A and/or CDKN2B  deletion
  - This is a secondary co-occurring alteration with prognostic significance. 


### Inclusion/exclusion criteria

Samples are _included_ for subtyping if we detect the following strings in the `pathology_diagnosis` field of `pbta-histologies.tsv`:

```
low-grade glioma/astrocytoma
ganglioglioma
```

Samples are _excluded_ if we detect the following strings in the `pathology_diagnosis` field of `pbta-histologies.tsv`:

```
dysembryoplastic neuroepithelial tumor
```

When `pathology_diagnosis == "Low-grade glioma/astrocytoma (WHO grade I/II)"`, we _exclude_ samples if we detect the following strings in `pathology_free_text_diagnosis`:

```
desmoplastic infantile astrocytoma
glioneuronal
```

These strings are stored in `lgat-subset/lgat_subtyping_path_dx_strings.json`, which is the output of the `00-LGAT-select-pathology-dx.R`. 

Prior to `release-v17-20200908`, we used `short_histology == "LGAT"` as our inclusion criteria.


### Preprocessing

The files in the `lgat-subset` were generated via `01-subset-files-for-LGAT.R` using files from the `release-v18-20201123`. If running locally using [the instructions below](#run-script), new `lgat-subset` files will be generated using the files symlinked in `data/`.

### Inputs from data download

* `pbta-histologies.tsv`: is used to subset samples according to [the criteria above](#inclusion-exclusion-criteria)
* `pbta-snv-consensus-mutation.maf.tsv.gz`: 

### Run script

```sh
bash run_subtyping.sh
```

This does not run the `00-v17-LGAT-select-pathology-dx` notebook, as that is intended to be run once and tied to a specific release (`release-v17-20200908`).

#### Order of scripts in subtyping

`01-subset-files-for-LGAT.Rmd`: generates subset of wgs LGAT sample annotated by subtypes snv mutation status as per following description:

columnname  | description | values
 --- | --- | ---
NF1_mut | somatic loss of NF1 via either missense, nonsense mutation | "Yes" mutation exists, "No" mutation is absent 
BRAF_V600E_mut | contains BRAF V600E or V599 SNV or non-canonical BRAF alterations such as p.V600ins or p.D594N | "Yes" mutation exists, "No" mutation is absent
MAPK_mut | contains mutation in KRAS, NRAS, HRAS, MAP2K1, MAP2K2, MAP2K1, ARAF, RAF1, BRAF (other than BRAF_V600E_mut) SNV or indel | "Yes" mutation exists, "No" mutation is absent
RTK_mut | harbors a MET,KIT or PDGFRA SNV | "Yes" mutation exists, "No" mutation is absent
FGFR_mut | harbors FGFR1 p.N546K, p.K656E, p.N577, or p. K687 hotspot mutations | "Yes" mutation exists, "No" mutation is absent
IDH_mut | harbors an IDH R132 mutation | "Yes" mutation exists, "No" mutation is absent
H3.3_mut | harbors an H3F3A K28M or G35R/V mutation | "Yes" mutation exists, "No" mutation is absent
H3.1_mut | harbors an HIST1H3B K28M|or HIST1H3C  K28M | "Yes" mutation exists, "No" mutation is absent

`02-subset-fusion-files-LGAT.Rmd`: generates subset of rna LGAT sample annotated by subtypes fusion status as per following description:

columnname | description | values
--- | --- | ---
KIAA_BRAF_fus | contains KIAA1549-BRAF fusion | "Yes" fusion exists, "No" fusion is absent
MAPK_fus | contains non-canonical BRAF fusion other than KIAA1549-BRAF; contains RAF1 fusion | "Yes" fusion exists, "No" fusion is absent
RTK_fus | harbors a fusion in ALK, ROS1, NTRK1, NTRK2,NTRK3 or PDGFR | "Yes" fusion exists, "No" fusion is absent
FGFR_fus | harbors FGFR1-TACC1 or other FGFR1 or FGFR2 fusion | "Yes" fusion exists, "No" fusion is absent
MYB_fus | harbors either a MYB-QKI fusion or other MYB or MYBL1 fusion | "Yes" fusion exists, "No" fusion is absent
 
`03-subset-cnv-files-LGAT.Rmd`: generates subset of dna LGAT sample annotated by subtypes defined by CNV as per the following description:

columnname | description | values
--- | --- | ---
FGFR_DUP_TANDEM | contains FGFR1 tandem duplication from SV | "Yes" fusion exists, "No" fusion is absent
FGFR_DUP | contains FGFR1 duplication event from consensus cnv seg | "Yes" fusion exists, "No" fusion is absent
CDKN2A_DEL | contains CDKN2A deletion | "Yes" fusion exists, "No" fusion is absent
CDKN2B_DEL | contains CDKN2B deletion | "Yes" fusion exists, "No" fusion is absent

