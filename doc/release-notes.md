# release notes
## curremt release
- release date: 2021-06-29
- status: available
- changes:
   - Within `histologies.tsv`:
     - Updated `cancer_group` logic to fix `Anaplastic` typo, remove benign tumors, and condense meningiomas [ticket48](https://github.com/PediatricOpenTargets/ticket-tracker/issues/48)
     - combine CBTN+PNOC into `cohort == PBTA` [ticket79](https://github.com/PediatricOpenTargets/ticket-tracker/issues/79)
     - GMKF tumor ploidy was added [ticket46](https://github.com/PediatricOpenTargets/ticket-tracker/issues/46)
     - harmonized tumor_descriptor per [ticket61](https://github.com/PediatricOpenTargets/ticket-tracker/issues/61) 
     - updated clinical info for NBL samples which were missing in source files per [ticket43](https://github.com/PediatricOpenTargets/ticket-tracker/issues/43)
     - updated `experimental_strategy` for targeted capture samples per [ticket62](https://github.com/PediatricOpenTargets/ticket-tracker/issues/62)
   - Add cnv files with PBTA+GMKF samples per [ticket44](https://github.com/PediatricOpenTargets/ticket-tracker/issues/44):
      - cnv-cnvkit.seg.gz
      - cnv-controlfreec.tsv.gz
      - cnv-consensus.seg.gz 
      - consensus_seg_annotated_cn_autosomes.tsv.gz
      - consensus_seg_annotated_cn_autosomes_xy.tsv.gz
    - Add EFO and MONDO cancer mapping file [ticket78](https://github.com/PediatricOpenTargets/ticket-tracker/issues/78):
    - Add ENSG to HUGO mapping file with RMTL designation [ticket84](https://github.com/PediatricOpenTargets/ticket-tracker/issues/84) and [ticket56](https://github.com/PediatricOpenTargets/ticket-tracker/issues/56)


## archived release
- release date: 2021-06-17
- status: available
- changes:
  - Removed
     - gtex_target_tcga-gene-expression-rsem-tpm-collapsed.polya.rds
     - gtex_target_tcga-gene-counts-rsem-expected_count-collapsed.rds 
  - Update RSEM files to include GTEXv8 files
     - gene-expression-rsem-tpm-collapsed.rds
     - gene-counts-rsem-expected_count-collapsed.rds
  - Update manifest to include GTEx manifest are now v8, TCGA manifest are from GDC and TARGET manifest from Diskin Lab and @afarrel
     - histologies.tsv
  - Added snv PBTA+GMKF maf file
     - snv-consensus-plus-hotspots.maf.tsv.gz 
  - Released files below

```
v5
├── WGS.hg38.lancet.300bp_padded.bed
├── WGS.hg38.lancet.unpadded.bed
├── WGS.hg38.mutect2.vardict.unpadded.bed
├── WGS.hg38.strelka2.unpadded.bed
├── WGS.hg38.vardict.100bp_padded.bed
├── fusion-arriba.tsv.gz
├── fusion-starfusion.tsv.gz
├── gene-counts-rsem-expected_count-collapsed.rds
├── gene-counts-rsem-expected_count.rds
├── gene-expression-rsem-tpm-collapsed.rds
├── gene-expression-rsem-tpm.rds
├── histologies.tsv
├── intersect_cds_lancet_strelka_mutect_WGS.bed
├── intersect_strelka_mutect_WGS.bed
├── snv-consensus-plus-hotspots.maf.tsv.gz
└── release-notes.md

```


## archived release
- release date: 2021-06-01
- status: available
- changes:
  - Update rnseq file names to be generic (currently includes PBTA and KFNBL)
    - `gene-counts-rsem-expected_count*`
    - `gene-expression-rsem-tpm*`
    - `fusion*`
  - Added mereged gtex_target_tcga expected_count file
    - `gtex_target_tcga-gene-counts-rsem-expected_count-collapsed.rds`
  - Combined PBTA, KFNBL, TARGET, TCGA, GTEX histologies into `histologies.tsv`  
  - Released files below

```
v4
├── WGS.hg38.lancet.300bp_padded.bed
├── WGS.hg38.lancet.unpadded.bed
├── WGS.hg38.mutect2.vardict.unpadded.bed
├── WGS.hg38.strelka2.unpadded.bed
├── WGS.hg38.vardict.100bp_padded.bed
├── fusion-arriba.tsv.gz
├── fusion-starfusion.tsv.gz
├── gene-counts-rsem-expected_count-collapsed.rds
├── gene-counts-rsem-expected_count.rds
├── gene-expression-rsem-tpm-collapsed.rds
├── gene-expression-rsem-tpm.rds
├── gtex_target_tcga-gene-counts-rsem-expected_count-collapsed.rds
├── histologies.tsv
├── intersect_cds_lancet_strelka_mutect_WGS.bed
├── intersect_strelka_mutect_WGS.bed
├── kfnbl-snv-lancet.vep.maf.gz
├── kfnbl-snv-mutect2.vep.maf.gz
├── kfnbl-snv-strelka2.vep.maf.gz
├── kfnbl-snv-vardict.vep.maf.gz
└── release-notes.md

```



## archived release
- release date: 2021-05-21
- status: available
- changes:
  - Update files due to addition of `SSF` column name 
    -  kfnbl-snv-vardict.vep.maf.gz
  - Added files for downstream analysis output files from `collapse-rnaseq`
    -  kfnbl-gene-counts-rsem-expected_count-collapsed.stranded.rds   
  - Released files below:

```
v3
├── WGS.hg38.lancet.300bp_padded.bed
├── WGS.hg38.lancet.unpadded.bed
├── WGS.hg38.mutect2.vardict.unpadded.bed
├── WGS.hg38.strelka2.unpadded.bed
├── WGS.hg38.vardict.100bp_padded.bed
├── gtex-gene-expression-rsem-tpm-collapsed.polya.rds
├── gtex-histologies.tsv
├── intersect_cds_lancet_strelka_mutect_WGS.bed
├── intersect_strelka_mutect_WGS.bed
├── kfnbl-fusion-arriba.tsv.gz
├── kfnbl-fusion-starfusion.tsv.gz
├── kfnbl-gene-counts-rsem-expected_count-collapsed.stranded.rds
├── kfnbl-gene-counts-rsem-expected_count.stranded.rds
├── kfnbl-gene-expression-kallisto.stranded.rds
├── kfnbl-gene-expression-rsem-fpkm-collapsed.stranded.rds
├── kfnbl-gene-expression-rsem-fpkm.stranded.rds
├── kfnbl-gene-expression-rsem-tpm-collapsed.stranded.rds
├── kfnbl-gene-expression-rsem-tpm.stranded.rds
├── kfnbl-histologies.tsv
├── kfnbl-isoform-counts-rsem-expected_count.stranded.rds
├── kfnbl-isoform-expression-rsem-tpm.stranded.rds
├── kfnbl-snv-lancet.vep.maf.gz
├── kfnbl-snv-mutect2.vep.maf.gz
├── kfnbl-snv-strelka2.vep.maf.gz
├── kfnbl-snv-vardict.vep.maf.gz
├── release-notes.md
├── target-gene-expression-rsem-tpm-collapsed.rds
├── target-histologies.tsv
├── tcga-gene-expression-rsem-tpm-collapsed.rds
└── tcga-histologies.tsv
```

## archived release
### release-v2
- release date: 2021-05-18
- status: available
- changes:
  - Updated rnaseq files with addition of BS_ETC8R0TD 
    - kfnbl-fusion-arriba.tsv.gz
    - kfnbl-fusion-starfusion.tsv.gz
    - kfnbl-gene-counts-rsem-expected_count.stranded.rds
    - kfnbl-gene-expression-kallisto.stranded.rds
    - kfnbl-gene-expression-rsem-fpkm.stranded.rds
    - kfnbl-gene-expression-rsem-fpkm-collapsed.stranded.rds
    - kfnbl-gene-expression-rsem-tpm.stranded.rds
    - kfnbl-gene-expression-rsem-tpm-collapsed.stranded.rds
    - kfnbl-isoform-counts-rsem-expected_count.stranded.rds
    - kfnbl-isoform-expression-rsem-tpm.stranded.rds
    - kfnbl-histologies.tsv 
  - Released files below: 
```
open-targets
└── v2
    ├── gtex-gene-expression-rsem-tpm-collapsed.polya.rds
    ├── gtex-histologies.tsv
    ├── intersect_cds_lancet_strelka_mutect_WGS.bed
    ├── intersect_strelka_mutect_WGS.bed
    ├── md5sum.txt
    ├── kfnbl-fusion-arriba.tsv.gz
    ├── kfnbl-fusion-starfusion.tsv.gz
    ├── kfnbl-gene-counts-rsem-expected_count.stranded.rds
    ├── kfnbl-gene-expression-kallisto.stranded.rds
    ├── kfnbl-gene-expression-rsem-fpkm.stranded.rds
    ├── kfnbl-gene-expression-rsem-fpkm-collapsed.stranded.rds
    ├── kfnbl-gene-expression-rsem-tpm.stranded.rds
    ├── kfnbl-gene-expression-rsem-tpm-collapsed.stranded.rds
    ├── kfnbl-isoform-counts-rsem-expected_count.stranded.rds
    ├── kfnbl-isoform-expression-rsem-tpm.stranded.rds
    ├── kfnbl-histologies.tsv    
    ├── kfnbl-snv-lancet.vep.maf.gz
    ├── kfnbl-snv-mutect2.vep.maf.gz
    ├── kfnbl-snv-strelka2.vep.maf.gz
    ├── kfnbl-snv-vardict.vep.maf.gz
    ├── release-notes.md 
    ├── target-gene-expression-rsem-tpm-collapsed.rds
    ├── target-histologies.tsv
    ├── tcga-gene-expression-rsem-tpm-collapsed.rds
    ├── tcga-histologies.tsv
    ├── WGS.hg38.lancet.300bp_padded.bed
    ├── WGS.hg38.lancet.unpadded.bed
    ├── WGS.hg38.mutect2.vardict.unpadded.bed
    ├── WGS.hg38.strelka2.unpadded.bed
    └── WGS.hg38.vardict.100bp_padded.bed
```

## archived release
### release-v1
- release date: 2021-04-21
- status: available
- changes:
  - Added files below: 
```
open-targets
└── v1
    ├── intersect_cds_lancet_strelka_mutect_WGS.bed
    ├── intersect_strelka_mutect_WGS.bed
    ├── md5sum.txt
    ├── kfnbl-gene-counts-rsem-expected_count.stranded.rds
    ├── kfnbl-gene-expression-rsem-fpkm.stranded.rds
    ├── kfnbl-gene-expression-rsem-tpm.stranded.rds
    ├── kfnbl-histologies.tsv    
    ├── kfnbl-snv-lancet.vep.maf.gz
    ├── kfnbl-snv-mutect2.vep.maf.gz
    ├── kfnbl-snv-strelka2.vep.maf.gz
    ├── kfnbl-snv-vardict.vep.maf.gz
    ├── release-notes.md
    ├── WGS.hg38.lancet.300bp_padded.bed
    ├── WGS.hg38.lancet.unpadded.bed
    ├── WGS.hg38.mutect2.vardict.unpadded.bed
    ├── WGS.hg38.strelka2.unpadded.bed
    └── WGS.hg38.vardict.100bp_padded.bed
```
