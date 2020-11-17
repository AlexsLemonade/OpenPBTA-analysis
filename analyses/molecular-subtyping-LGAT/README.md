# Molecular Subtyping of Low Grade Gliomaa

**Module authors:** Krutika Gaonkar ([@kgaonkar6](https://github.com/kgaonkar6)), Jaclyn Taroni ([@jaclyn-taroni](https://github.com/jaclyn-taroni))

In this analysis we subtype LGAT samples according to the presence/absence of _BRAF_ fusions and BRAF V600E point mutations. 

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

These strings are stored in `lgat-subset/lgat_subtyping_path_dx_strings.json`, which is the output of the `00-LGAT-select-pathology-dx` notebook. 

Prior to `release-v17-20200908`, we used `short_histology == "LGAT"` as our inclusion criteria.


### Preprocessing

The files in the `lgat-subset` were generated via `01-subset-files-for-LGAT.R` using files from the `release-v17-20200908`. If running locally using [the instructions below](#run-script), new `lgat-subset` files will be generated using the files symlinked in `data/`.

### Inputs from data download

* `pbta-histologies.tsv`: is used to subset samples according to [the criteria above](#inclusion-exclusion-criteria)
* `pbta-snv-consensus-mutation.maf.tsv.gz` 
* `fusion_summary_lgat_foi.tsv` subset of fusion as per LGAT subtype  

### Run script

```sh
bash run_subtyping.sh
```

This does not run the `00-LGAT-select-pathology-dx` notebook, as that is intended to be run once and tied to a specific release (`release-v17-20200908`).

#### Order of scripts in subtyping

`01-subset-files-for-LGAT.R`: generates subset of wgs LGAT sample annotated by subtypes snv mutation status as per following description:

columnname  | description | values
 --- | --- | ---
NF1_mut | somatic loss of NF1 via either missense, nonsense mutation | "Yes" mutation exists, "No" mutation is absent 
BRAF_V600E_mut | contains BRAF V600E or V599 SNV or non-canonical BRAF alterations such as p.V600ins or p.D594N | "Yes" mutation exists, "No" mutation is absent
MAPK_mut | contains mutation in KRAS, NRAS, HRAS, MAP2K1, MAP2K2, MAP2K1, ARAF SNV or indel | "Yes" mutation exists, "No" mutation is absent
RTK_mut | harbors a MET,KIT or PDGFRA SNV | "Yes" mutation exists, "No" mutation is absent
FGFR_mut | harbors FGFR1 p.N546K, p.K656E, p.N577, or p. K687 hotspot mutations | "Yes" mutation exists, "No" mutation is absent
IDH_mut | harbors an IDH R132 mutation | "Yes" mutation exists, "No" mutation is absent
H3F3A_mut | harbors an H3F3A K28M or G35R/V mutation | "Yes" mutation exists, "No" mutation is absent
HIST1H3B_mut | harbors an HIST1H3B K28M | "Yes" mutation exists, "No" mutation is absent
HIST1H3C_mut | harbors and HIST1H3C  K28M | "Yes" mutation exists, "No" mutation is absent
 
`02-subset-fusion-files-LGAT.R`: generates subset of rna LGAT sample annotated by subtypes fusion status as per following description:
columnname | description | values
--- | --- | ---
KIAA_BRAF_fus | contains KIAA1549-BRAF fusion | "Yes" fusion exists, "No" fusion is absent
MAPK_fus | contains non-canonical BRAF fusion other than KIAA1549-BRAF; contains RAF1 fusion | "Yes" fusion exists, "No" fusion is absent
RTK_fus | harbors a fusion in ALK, ROS1, NTRK1, NTRK2,NTRK3 or PDGFR | "Yes" fusion exists, "No" fusion is absent
FGFR_fus | harbors FGFR1-TACC1 or other FGFR1 or FGFR2 fusion | "Yes" fusion exists, "No" fusion is absent
MYB_fus | harbors either a MYB-QKI fusion or other MYB or MYBL1 fusion | "Yes" fusion exists, "No" fusion is absent
