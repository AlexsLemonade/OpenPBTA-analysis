# release notes

## current release
- release date: 2021-09-01
- status: available
- overview of changes:
  - This particular release is just an update of `efo-mondo-map.tsv` per discussion in [ticket 182](https://github.com/PediatricOpenTargets/ticket-tracker/issues/182).
  - MONDO and EFO codes were manually reviewed and assigned to `cancer_group` in OpenPedCan.  

```
v9
├── WGS.hg38.lancet.300bp_padded.bed
├── WGS.hg38.lancet.unpadded.bed
├── WGS.hg38.mutect2.vardict.unpadded.bed
├── WGS.hg38.strelka2.unpadded.bed
├── WGS.hg38.vardict.100bp_padded.bed
├── cnv-cnvkit.seg.gz
├── cnv-controlfreec.tsv.gz
├── consensus_wgs_plus_cnvkit_wxs_autosomes.tsv.gz
├── consensus_wgs_plus_cnvkit_wxs_x_and_y.tsv.gz
├── consensus_wgs_plus_cnvkit_wxs.tsv.gz
├── data-files-description.md
├── efo-mondo-map.tsv
├── ensg-hugo-rmtl-mapping.tsv
├── fusion-arriba.tsv.gz
├── fusion-starfusion.tsv.gz
├── fusion-putative-oncogenic.tsv
├── gene-counts-rsem-expected_count-collapsed.rds
├── gene-expression-rsem-tpm-collapsed.rds
├── histologies.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.tsv
├── independent-specimens.rnaseq.primary.eachcohort.tsv
├── independent-specimens.rnaseq.relapse.eachcohort.tsv
├── independent-specimens.wgswxspanel.primary.tsv
├── independent-specimens.wgswxspanel.relapse.tsv
├── independent-specimens.rnaseq.primary.tsv
├── independent-specimens.rnaseq.relapse.tsv
├── intersect_cds_lancet_strelka_mutect_WGS.bed
├── intersect_strelka_mutect_WGS.bed
├── release-notes.md
├── snv-consensus-plus-hotspots.maf.tsv.gz
├── sv-manta.tsv.gz
├── uberon-map-gtex-group.tsv
└── uberon-map-gtex-subgroup.tsv

```

## archived release
- release date: 2021-08-20
- status: available
- overview of changes:
    - This particular release is mainly to include 412 tumor/normal pairs of TARGET WXS samples as listed in [ticket 111](https://github.com/PediatricOpenTargets/ticket-tracker/issues/111). Detailed changes see below. 
- detailed changes: 
    - Histology file updates - the master ticket is d3b center [ticket 43](https://github.com/d3b-center/D3b-codes/pull/43):
      - 412 tumor/normal TARGET WXS samples were included in the histologies.tsv file per [ticket 111](https://github.com/PediatricOpenTargets/ticket-tracker/issues/111)v
      - `sample_id` and `aliquot_id` for some TARGET samples were miscoded before. This release fixed the issue per [ticket 145](https://github.com/PediatricOpenTargets/ticket-tracker/issues/145)
      - `sample_id` and `aliquot_id` were updated using the GTEx coding nomenclature GTEX-[donor ID]-[tissue site ID]-SM-[aliquot ID] (https://www.gtexportal.org/home/faq#sampleIdFormat)
      - `primary_site` for GTEx samples were updated to match `gtex_subgroup` column, changing `Whole Blood` to `Blood` and `Brain - Cerebellar Hemisphere` to `Brain - Cerebellum`
      - `broad_histology` for TARGET samples were updated as following: `Acute Lymphoblastic Leukemia` and `Acute Myeloid Leukemia` were merged to `Hematologic malignancy`; `Clear cell sarcoma of the kidney`, `Rhabdoid tumor`, and `Wilms tumor` were combined as `Renal tumor`; `Osteosarcoma` is changed to `Mesenchymal non-meningothelial tumor` and `Neuroblastoma` is converted to `Embryonal tumor`. See [ticket 136](https://github.com/PediatricOpenTargets/ticket-tracker/issues/136) and [ticket 176] (https://github.com/PediatricOpenTargets/ticket-tracker/issues/176)
      - For `gtex_group == "Cells"`, the `composition` column is changed from `Solid Tissue` to `Derived Cell Line` per discussion in d3b center [ticket 43](https://github.com/d3b-center/D3b-codes/pull/43)
      - For `short_histology`,  neuroblastoma samples previously annotated as `NBL` or `Embryonal tumor` are converted to `Neuroblastoma` to be consistent with other samples 
      - Updated MB subtypes in the histologies file per [ticket 148](https://github.com/PediatricOpenTargets/ticket-tracker/issues/148)
      - Some GMKF WGS samples does not have `germline_sex_estimate` as indicated in [ticket 168](https://github.com/PediatricOpenTargets/ticket-tracker/issues/168). Added using the file in the discussion of [ticket 159](https://github.com/PediatricOpenTargets/ticket-tracker/issues/159)
      - Some TARGET WXS samples miss `tumor_ploidy` that are actually available - and those samples now have updated ploidy using the file in the discussion session of [ticket 160](https://github.com/PediatricOpenTargets/ticket-tracker/issues/160)
    
    - Update DNA related files to include TARGET WXS DNA (412 tumor/normal pairs):
      - cnv-cnvkit.seg.gz [ticket 156](https://github.com/PediatricOpenTargets/ticket-tracker/issues/156)
      - cnv-controlfreec.tsv.gz [ticket 156](https://github.com/PediatricOpenTargets/ticket-tracker/issues/156)
      - snv-consensus-plus-hotspots.maf.tsv.gz [ticket 156](https://github.com/PediatricOpenTargets/ticket-tracker/issues/156)
    - Update method to call CNV consensus (WGS) as described in [ticket 134](https://github.com/PediatricOpenTargets/ticket-tracker/issues/134) and [ticket 149](https://github.com/PediatricOpenTargets/ticket-tracker/issues/149). Briefly, CNV called by `MantaSV` were filtered to contain only `filter == "PASS"` before going into the consensus calling workflow. In the subsequent step, instead of only retaining CNV calls that have 50% reciprocal overlap between callers (which was too stringent), the criteria is expanded to include small CNV regions that are 90% covered by a larger CNV. The consensus is the overlapping region. 
      - cnv_consensus_seg.gz (WGS samples only - in S3 bucket s3://kf-openaccess-us-east-1-prd-pbta/open-targets/v8/ but not in md5sum.txt file)
    - As a result of changing consensus calling criteria and adding new TARGET WXS DNA sample results to `cnv-cnvkit.seg.gz` and `cnv-controlfreec.tsv.gz`, the following files were updated:
      - consensus_wgs_plus_cnvkit_wxs_autosomes.tsv.gz [ticket 159](https://github.com/PediatricOpenTargets/ticket-tracker/issues/159) and [ticket 160](https://github.com/PediatricOpenTargets/ticket-tracker/issues/160) (only in S3 bucket s3://kf-openaccess-us-east-1-prd-pbta/open-targets/v8/- not included in automatic download)
      - consensus_wgs_plus_cnvkit_wxs_x_and_y.tsv.gz [ticket 159](https://github.com/PediatricOpenTargets/ticket-tracker/issues/159) and [ticket 160](https://github.com/PediatricOpenTargets/ticket-tracker/issues/160) (only in S3 bucket s3://kf-openaccess-us-east-1-prd-pbta/open-targets/v8/- not included in automatic download)
      
    - Added `consensus_wgs_plus_cnvkit_wxs.tsv.gz` which is a merge of `consensus_wgs_plus_cnvkit_wxs_autosomes.tsv.gz` and `consensus_wgs_plus_cnvkit_wxs_x_and_y.tsv.gz` per [ticket 161](https://github.com/PediatricOpenTargets/ticket-tracker/issues/161)
    - Updated `ensg-hugo-rmtl-mapping.tsv` file per [ticket 146](https://github.com/PediatricOpenTargets/ticket-tracker/issues/146). The previous release of this file does not contain all gene ENSG IDs and symbols that are present in `snv-consensus-plus-hotspots.maf.tsv.gz`. This update merged GENCODE V28 and V38 to allow inclusion of more gene ENSG IDs and symbols. 
    - Futher update `ensg-hugo-rmtl-mapping.tsv` [PR 48 D3b codes](https://github.com/d3b-center/D3b-codes/pull/48) to include all gene ENSG ID to symbol mappings in v7 `ensg-hugo-rmtl-mapping.tsv`.

    - Update independent samples files to include TARGET WXS DNA (412 tumor/normal pairs) - [ticket 165](https://github.com/PediatricOpenTargets/ticket-tracker/issues/165). Previously, these files were not added to our releases. Starting this release, we will also add independent sample list to our release as well. 
    - Updated independent samples so that the `Kids_First_Biospecimen_ID` for `allcohorts` and `eachcohort` match if possible: [ticket 135](https://github.com/PediatricOpenTargets/ticket-tracker/issues/135)
    - Updated independent sample module to arrange by `Kids_First_Biospecimen_ID` before writing out the file: [ticket 179](https://github.com/PediatricOpenTargets/ticket-tracker/issues/179)
    - For now, we will add files that are used by analyses modules in this file and these are the following: 
      - independent-specimens.wgswxspanel.primary.eachcohort.tsv
      - independent-specimens.wgswxspanel.relapse.eachcohort.tsv
      - independent-specimens.rnaseq.primary.eachcohort.tsv
      - independent-specimens.rnaseq.relapse.eachcohort.tsv
      - independent-specimens.wgswxspanel.primary.tsv
      - independent-specimens.wgswxspanel.relapse.tsv
      - independent-specimens.rnaseq.primary.tsv
      - independent-specimens.rnaseq.relapse.tsv
    - Updated `fusion-putative-oncogenic.tsv` since at the last step of putative oncogenic fusion filtering, we filter out fusions seen in > 4 broad_histology since they are likely artifacts and with the update of broad_histology, the result will be updated [ticket 175](https://github.com/PediatricOpenTargets/ticket-tracker/issues/175)
      - fusion-putative-oncogenic.tsv

```
v8
├── WGS.hg38.lancet.300bp_padded.bed
├── WGS.hg38.lancet.unpadded.bed
├── WGS.hg38.mutect2.vardict.unpadded.bed
├── WGS.hg38.strelka2.unpadded.bed
├── WGS.hg38.vardict.100bp_padded.bed
├── cnv-cnvkit.seg.gz
├── cnv-controlfreec.tsv.gz
├── consensus_wgs_plus_cnvkit_wxs_autosomes.tsv.gz
├── consensus_wgs_plus_cnvkit_wxs_x_and_y.tsv.gz
├── consensus_wgs_plus_cnvkit_wxs.tsv.gz
├── data-files-description.md
├── efo-mondo-map.tsv
├── ensg-hugo-rmtl-mapping.tsv
├── fusion-arriba.tsv.gz
├── fusion-starfusion.tsv.gz
├── fusion-putative-oncogenic.tsv
├── gene-counts-rsem-expected_count-collapsed.rds
├── gene-expression-rsem-tpm-collapsed.rds
├── histologies.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.tsv
├── independent-specimens.rnaseq.primary.eachcohort.tsv
├── independent-specimens.rnaseq.relapse.eachcohort.tsv
├── independent-specimens.wgswxspanel.primary.tsv
├── independent-specimens.wgswxspanel.relapse.tsv
├── independent-specimens.rnaseq.primary.tsv
├── independent-specimens.rnaseq.relapse.tsv
├── intersect_cds_lancet_strelka_mutect_WGS.bed
├── intersect_strelka_mutect_WGS.bed
├── release-notes.md
├── snv-consensus-plus-hotspots.maf.tsv.gz
├── sv-manta.tsv.gz
├── uberon-map-gtex-group.tsv
└── uberon-map-gtex-subgroup.tsv

```

- release date: 2021-07-23
- status: available
- changes:
   - Updated EFO/MONDO mapping file per [ticket 88](https://github.com/PediatricOpenTargets/ticket-tracker/issues/88)
   - Added GTEX UBERON mapping files for subgroup and group per [ticket 85](https://github.com/PediatricOpenTargets/ticket-tracker/issues/85)
     - Collapsed `Cerebellum hemisphere` and `Cerebellum` to `Cerebellum` since GTEX has the same UBERON code listed for both per [ticket 106](https://github.com/PediatricOpenTargets/ticket-tracker/issues/106)
     - Renamed `Whole Blood` to `Blood` per John Maris's suggestion to alphabetize and [ticket 106](https://github.com/PediatricOpenTargets/ticket-tracker/issues/106). 
       - Note: v8 gtex lists "whole blood" subgroup and "whole blood" group as UBERON_0013756, which is mapped to venous blood and "Whole blood" maps to UBERON_0000178. 
       - After inquiry at GTEx, we were told they are equivalent terms as seen in [this link](https://www.ebi.ac.uk/ols/ontologies/uberon/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FUBERON_0013756)
    - Updated `ensg-hugo-rmtl-v1-mapping.tsv` with minor updates according to [ticket 125](https://github.com/PediatricOpenTargets/ticket-tracker/issues/125) and changed filename to `ensg-hugo-rmtl-mapping.tsv`
    - Histology file updates:
      - Collapsed `Cerebellum hemisphere` and `Cerebellum` to `Cerebellum` since GTEX has the same UBERON code listed for both per [ticket 106](https://github.com/PediatricOpenTargets/ticket-tracker/issues/106)
      - Renamed `Whole Blood` to `Blood` per John Maris's suggestion to alphabetize and [ticket 106](https://github.com/PediatricOpenTargets/ticket-tracker/issues/106). 
      - Added inferred strandedness to RNA-Seq samples per [ticket 104](https://github.com/PediatricOpenTargets/ticket-tracker/issues/104)
      - Added `broad_tumor_descriptor` to designate grouped `Diagnosis` and `Relapse` samples used in SNV, CNV, Fusion tables as well as for grouping on pedcbio per [ticket 109](https://github.com/PediatricOpenTargets/ticket-tracker/issues/109)
      - Added ploidy information for TARGET AML and NBL WXS samples per [ticket 121](https://github.com/PediatricOpenTargets/ticket-tracker/issues/121)
      - Collapsed TARGET ids containing suffixes to match the BAM file sample IDs from GDC and match the RDS processed files per [comment here](https://github.com/d3b-center/D3b-codes/pull/41#issuecomment-885809293)
    - Added TARGET NBL and AML WXS, PBTA WXS CNV calls to `cnv-cnvkit.seg.gz` and `cnv-controlfreec.tsv.gz` per [ticket 80](https://github.com/PediatricOpenTargets/ticket-tracker/issues/80) 
    - Added `consensus_wgs_plus_cnvkit_wxs_x_and_y.tsv.gz` (removed WGS only file `consensus_seg_annotated_cn_autosomes.tsv.gz`) and `consensus_wgs_plus_cnvkit_wxs_autosomes.tsv.gz` (removed WGS only file `consensus_seg_annotated_cn_x_and_y.tsv.gz`) containing consensus WGS and CNVkit WXS data per [ticket 102](https://github.com/PediatricOpenTargets/ticket-tracker/issues/102)
    - Updated RNA-Seq files to include TARGET RNA (N = 1329) samples:
      - fusion-arriba.tsv.gz
      - fusion-starfusion.tsv.gz
      - fusion-putative-oncogenic.tsv
      - gene-counts-rsem-expected_count-collapsed.rds
      - gene-expression-rsem-tpm-collapsed.rds

## archived release
- release date: 2021-06-29
- status: available
- changes:
   - Within `histologies.tsv`:
     - Updated `cancer_group` logic to make updates as per [ticket48](https://github.com/PediatricOpenTargets/ticket-tracker/issues/48)
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
