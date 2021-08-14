# release notes
## current release
### release-v20-20210726
- release data: 2021-07-26
- status: availale
- changes:
 - Update copy_number_consensus_call with manta calls with FILTER==PASS [#1113](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/1113)
 - Update reciprocal overlap rules for consensus CNV [#1125](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/1125)
    - pbta-cnv-consensus.seg.gz
    - pbta-cnv-consensus-gistic.zip
 - Update focal-cn-file-preparation with deep deletion status and with updated pbta-cnv-consensus.seg.gz above [#1099](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/1099)
    - consensus_seg_annotated_cn_autosomes.tsv.gz
    - consensus_seg_annotated_cn_x_and_y.tsv.gz
 - Add medulloblastoma subtypes being extracted from pathology reports [#746](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/746)
    - pbta-mb-pathology-subtypes.tsv
 - Add cnv and indel hostspots from MSKCC and TERT promoters [#932](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/932) 
    - pbta-snv-scavenged-hotspots.maf.tsv.gz
 - Update histology information 
   - latest pathology diagnosis for 7316-2686 [#1032](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/1032)
   - BS_0TCRV9AC as sample_type == Normal
   - Add `cancer_group` for plot labeling per [#917](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/917)
- folder structure
```
data
└── release-v19-20210423
    ├── consensus_seg_annotated_cn_autosomes.tsv.gz
    ├── consensus_seg_annotated_cn_x_and_y.tsv.gz
    ├── data-files-description.md
    ├── fusion_summary_embryonal_foi.tsv
    ├── fusion_summary_ependymoma_foi.tsv
    ├── fusion_summary_ewings_foi.tsv
    ├── fusion_summary_lgat_foi.tsv
    ├── independent-specimens.rnaseq.primary-plus-polya.tsv
    ├── independent-specimens.rnaseq.primary-plus-stranded.tsv
    ├── independent-specimens.wgs.primary-plus.tsv
    ├── independent-specimens.wgs.primary.tsv
    ├── independent-specimens.wgswxs.primary-plus.tsv
    ├── independent-specimens.wgswxs.primary.tsv
    ├── intersect_cds_lancet_strelka_mutect_WGS.bed
    ├── intersect_cds_lancet_WXS.bed
    ├── intersect_strelka_mutect_WGS.bed
    ├── md5sum.txt
    ├── pbta-cnv-cnvkit-gistic.zip
    ├── pbta-cnv-cnvkit.seg.gz
    ├── pbta-cnv-consensus-gistic.zip
    ├── pbta-cnv-consensus.seg.gz
    ├── pbta-cnv-controlfreec.tsv.gz
    ├── pbta-fusion-arriba.tsv.gz
    ├── pbta-fusion-putative-oncogenic.tsv
    ├── pbta-fusion-recurrently-fused-genes-byhistology.tsv
    ├── pbta-fusion-recurrently-fused-genes-bysample.tsv
    ├── pbta-fusion-starfusion.tsv.gz
    ├── pbta-gene-counts-rsem-expected_count.polya.rds
    ├── pbta-gene-counts-rsem-expected_count.stranded.rds
    ├── pbta-gene-expression-kallisto.polya.rds
    ├── pbta-gene-expression-kallisto.stranded.rds
    ├── pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
    ├── pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
    ├── pbta-gene-expression-rsem-fpkm.polya.rds
    ├── pbta-gene-expression-rsem-fpkm.stranded.rds
    ├── pbta-gene-expression-rsem-tpm.polya.rds
    ├── pbta-gene-expression-rsem-tpm.stranded.rds
    ├── pbta-histologies.tsv    
    ├── pbta-histologies-base.tsv
    ├── pbta-isoform-counts-rsem-expected_count.polya.rds
    ├── pbta-isoform-counts-rsem-expected_count.stranded.rds
    ├── pbta-isoform-expression-rsem-tpm.polya.rds
    ├── pbta-isoform-expression-rsem-tpm.stranded.rds
    ├── pbta-mend-qc-manifest.tsv
    ├── pbta-mend-qc-results.tar.gz
    ├── pbta-snv-consensus-mutation.maf.tsv.gz
    ├── pbta-snv-consensus-hotspot-mutation.maf.tsv.gz
    ├── pbta-snv-consensus-mutation-tmb-all.tsv
    ├── pbta-snv-consensus-mutation-tmb-coding.tsv
    ├── pbta-snv-lancet.vep.maf.gz
    ├── pbta-snv-mutect2.vep.maf.gz
    ├── pbta-snv-scavenged-hotspots.maf.tsv.gz
    ├── pbta-snv-strelka2.vep.maf.gz
    ├── pbta-snv-vardict.vep.maf.gz
    ├── pbta-star-log-final.tar.gz
    ├── pbta-star-log-manifest.tsv
    ├── pbta-sv-manta.tsv.gz
    ├── pbta-tcga-manifest.tsv
    ├── pbta-tcga-snv-lancet.vep.maf.gz
    ├── pbta-tcga-snv-mutect2.vep.maf.gz
    ├── pbta-tcga-snv-strelka2.vep.maf.gz
    ├── release-notes.md
    ├── StrexomeLite_hg38_liftover_100bp_padded.bed
    ├── StrexomeLite_Targets_CrossMap_hg38_filtered_chr_prefixed.bed
    ├── WGS.hg38.lancet.300bp_padded.bed
    ├── WGS.hg38.lancet.unpadded.bed
    ├── WGS.hg38.mutect2.vardict.unpadded.bed
    ├── WGS.hg38.strelka2.unpadded.bed
    ├── WGS.hg38.vardict.100bp_padded.bed
    ├── WXS.hg38.100bp_padded.bed
    ├── WXS.hg38.lancet.400bp_padded.bed
    ├── intersected_whole_exome_agilent_designed_120_AND_tcga_6k_genes.Gh38.bed
    ├── intersected_whole_exome_agilent_plus_tcga_6k_AND_tcga_6k_genes.Gh38.bed
    ├── tcga_6k_genes.targetIntervals.Gh38.bed
    ├── tcga_6k_genes.targetIntervals.bed
    ├── tcga-snv-consensus-snv.maf.tsv.gz
    ├── tcga-snv-mutation-tmb-all.tsv
    ├── tcga-snv-mutation-tmb-coding.tsv
    ├── whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.Gh38.bed
    ├── whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.bed
    ├── whole_exome_agilent_designed_120.targetIntervals.Gh38.bed
    ├── whole_exome_agilent_designed_120.targetIntervals.bed
    ├── whole_exome_agilent_plus_tcga_6k.targetIntervals.Gh38.bed
    └── whole_exome_agilent_plus_tcga_6k.targetIntervals.bed
``` 


## archived release
### release-v19-20210423
- release date: 2021-04-23
- status: available
- changes:
  - Remove per BS_JXF8A2A6 due to mis-identification per [#862](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/862)
    - Updates:
      - pbta-fusion-starfusion.tsv.gz
      - pbta-fusion-arriba.tsv.gz
      - pbta-gene-expression-rsem-fpkm.stranded.rds
      - pbta-gene-counts-rsem-expected_count.stranded.rds
      - pbta-gene-expression-rsem-tpm.stranded.rds
      - pbta-gene-expression-kallisto.stranded.rds
      - pbta-isoform-counts-rsem-expected_count.stranded.rds
      - pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
      - pbta-isoform-expression-rsem-tpm.stranded.rds
      - independent-specimens.rnaseq.primary-plus-stranded.tsv
      - pbta-star-log-final.tar.gz
      - pbta-star-log-manifest.tsv
      - pbta-mend-qc-manifest.tsv
      - pbta-mend-qc-results.tar.gz
      - pbta-fusion-putative-oncogenic.tsv
      - pbta-fusion-recurrently-fused-genes-byhistology.tsv
      - pbta-fusion-recurrently-fused-genes-bysample.tsv
      - independent-specimens.rnaseq.primary-plus-stranded.tsv
      - fusion_summary_ewings_foi.tsv
      - fusion_summary_lgat_foi.tsv
      - fusion_summary_ependymoma_foi.tsv 
      - fusion_summary_embryonal_foi.tsv
  - Add INDEL renormalization of MAF files per [#1024](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/1024)
    - Updates:
      - pbta-snv-strelka2.vep.maf.gz
      - pbta-snv-lancet.vep.maf.gz            
      - pbta-snv-vardict.vep.maf.gz            
      - pbta-snv-mutect2.vep.maf.gz 
      - pbta-tcga-snv-strelka2.vep.maf.gz    
      - pbta-tcga-snv-mutect2.vep.maf.gz       
      - pbta-tcga-snv-lancet.vep.maf.gz 
  - Update pbta-histologies.tsv:
    - Add `PFS_days` to `pbta-histologies.tsv` file per [#963](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/963) 
    - Pull latest clinical data
    - Rerun molecular subtyping modules to date
    - Update PNOC003 `cohort` to "PNOC" for inclusion of later PNOC trial samples
    - Update `CNS_region` from manual review per [#1025](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/1025)
  - Add `pbta-mb-pathology-subtypes.tsv` containing pathology report molecular subtypes for medulloblastoma to be used as a "truth" set for medulloblastoma classifier testing per [#746](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/746)
- folder structure:
```
data
└── release-v19-20210423
    ├── consensus_seg_annotated_cn_autosomes.tsv.gz
    ├── consensus_seg_annotated_cn_x_and_y.tsv.gz
    ├── data-files-description.md
    ├── fusion_summary_embryonal_foi.tsv
    ├── fusion_summary_ependymoma_foi.tsv
    ├── fusion_summary_ewings_foi.tsv
    ├── fusion_summary_lgat_foi.tsv
    ├── independent-specimens.rnaseq.primary-plus-polya.tsv
    ├── independent-specimens.rnaseq.primary-plus-stranded.tsv
    ├── independent-specimens.wgs.primary-plus.tsv
    ├── independent-specimens.wgs.primary.tsv
    ├── independent-specimens.wgswxs.primary-plus.tsv
    ├── independent-specimens.wgswxs.primary.tsv
    ├── intersect_cds_lancet_strelka_mutect_WGS.bed
    ├── intersect_cds_lancet_WXS.bed
    ├── intersect_strelka_mutect_WGS.bed
    ├── md5sum.txt
    ├── pbta-cnv-cnvkit-gistic.zip
    ├── pbta-cnv-cnvkit.seg.gz
    ├── pbta-cnv-consensus-gistic.zip
    ├── pbta-cnv-consensus.seg.gz
    ├── pbta-cnv-controlfreec.tsv.gz
    ├── pbta-fusion-arriba.tsv.gz
    ├── pbta-fusion-putative-oncogenic.tsv
    ├── pbta-fusion-recurrently-fused-genes-byhistology.tsv
    ├── pbta-fusion-recurrently-fused-genes-bysample.tsv
    ├── pbta-fusion-starfusion.tsv.gz
    ├── pbta-gene-counts-rsem-expected_count.polya.rds
    ├── pbta-gene-counts-rsem-expected_count.stranded.rds
    ├── pbta-gene-expression-kallisto.polya.rds
    ├── pbta-gene-expression-kallisto.stranded.rds
    ├── pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
    ├── pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
    ├── pbta-gene-expression-rsem-fpkm.polya.rds
    ├── pbta-gene-expression-rsem-fpkm.stranded.rds
    ├── pbta-gene-expression-rsem-tpm.polya.rds
    ├── pbta-gene-expression-rsem-tpm.stranded.rds
    ├── pbta-histologies.tsv    
    ├── pbta-histologies-base.tsv
    ├── pbta-isoform-counts-rsem-expected_count.polya.rds
    ├── pbta-isoform-counts-rsem-expected_count.stranded.rds
    ├── pbta-isoform-expression-rsem-tpm.polya.rds
    ├── pbta-isoform-expression-rsem-tpm.stranded.rds
    ├── pbta-mend-qc-manifest.tsv
    ├── pbta-mend-qc-results.tar.gz
    ├── pbta-snv-consensus-mutation.maf.tsv.gz
    ├── pbta-snv-consensus-mutation-tmb-all.tsv
    ├── pbta-snv-consensus-mutation-tmb-coding.tsv
    ├── pbta-snv-lancet.vep.maf.gz
    ├── pbta-snv-mutect2.vep.maf.gz
    ├── pbta-snv-scavenged-hotspots.maf.tsv.gz
    ├── pbta-snv-strelka2.vep.maf.gz
    ├── pbta-snv-vardict.vep.maf.gz
    ├── pbta-star-log-final.tar.gz
    ├── pbta-star-log-manifest.tsv
    ├── pbta-sv-manta.tsv.gz
    ├── pbta-tcga-manifest.tsv
    ├── pbta-tcga-snv-lancet.vep.maf.gz
    ├── pbta-tcga-snv-mutect2.vep.maf.gz
    ├── pbta-tcga-snv-strelka2.vep.maf.gz
    ├── release-notes.md
    ├── StrexomeLite_hg38_liftover_100bp_padded.bed
    ├── StrexomeLite_Targets_CrossMap_hg38_filtered_chr_prefixed.bed
    ├── WGS.hg38.lancet.300bp_padded.bed
    ├── WGS.hg38.lancet.unpadded.bed
    ├── WGS.hg38.mutect2.vardict.unpadded.bed
    ├── WGS.hg38.strelka2.unpadded.bed
    ├── WGS.hg38.vardict.100bp_padded.bed
    ├── WXS.hg38.100bp_padded.bed
    ├── WXS.hg38.lancet.400bp_padded.bed
    ├── intersected_whole_exome_agilent_designed_120_AND_tcga_6k_genes.Gh38.bed
    ├── intersected_whole_exome_agilent_plus_tcga_6k_AND_tcga_6k_genes.Gh38.bed
    ├── tcga_6k_genes.targetIntervals.Gh38.bed
    ├── tcga_6k_genes.targetIntervals.bed
    ├── tcga-snv-consensus-snv.maf.tsv.gz
    ├── tcga-snv-mutation-tmb-all.tsv
    ├── tcga-snv-mutation-tmb-coding.tsv
    ├── whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.Gh38.bed
    ├── whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.bed
    ├── whole_exome_agilent_designed_120.targetIntervals.Gh38.bed
    ├── whole_exome_agilent_designed_120.targetIntervals.bed
    ├── whole_exome_agilent_plus_tcga_6k.targetIntervals.Gh38.bed
    └── whole_exome_agilent_plus_tcga_6k.targetIntervals.bed
```

## archived release
### release-v18-20201123
- release date: 2020-11-23
- status: available
- changes:
  - Add `pbta-histologies-base.tsv` file, which does not contain `molecular_subtype` nor `integrated_diagnosis` and will be used for molecular subtyping modules in conjunction with data releases
  - Add RNA independent specimen lists from [#795](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/795) and [#797](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/797)
  - Update DNA independent specimen lists per comment [here](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/750#issuecomment-691526896)
  - Add LGAT fusion summary file from [#808](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/808) and [#830](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/830)
  - Update fusion files per comment [here](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/750#issuecomment-691529675)
- Add kinase domain and reciprocal information to pbta-fusion-putative-oncogenic.tsv per [#812](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/812) and [#816](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/816) and [#821](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/821)
  - Add 8 rRNA-depleted, total stranded RNA-Seq for tumors which previously had polyA sequencing performed: BS_FXJY0MNH, BS_KABQQA0T, BS_D7XRFE0R, BS_SHJA4MR0, BS_HE0WJRW6, BS_8QB4S4VA, BS_FN07P04C, BS_SB12W1XT per [#749](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/749)
  - Update pbta-histologies.tsv:
    - Add new RNA samples above
    - Pull latest clinical data
    - Add extent of tumor resection per [comment](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/750#issue-688364970)
    - Rerun molecular subtyping modules to date
    - Harmonize PNOC003 `pathology_diagnosis` to CBTN `pathology_diagnosis` and move original diagnosis to `pathology_free_text_diagnosis`
    - Update CBTTC `cohort` to "CBTN"
- folder structure:
```
data
└── release-v18-20201123
    ├── consensus_seg_annotated_cn_autosomes.tsv.gz
    ├── consensus_seg_annotated_cn_x_and_y.tsv.gz
    ├── data-files-description.md
    ├── fusion_summary_embryonal_foi.tsv
    ├── fusion_summary_ependymoma_foi.tsv
    ├── fusion_summary_ewings_foi.tsv
    ├── fusion_summary_lgat_foi.tsv
    ├── independent-specimens.rnaseq.primary-plus-polya.tsv
    ├── independent-specimens.rnaseq.primary-plus-stranded.tsv
    ├── independent-specimens.wgs.primary-plus.tsv
    ├── independent-specimens.wgs.primary.tsv
    ├── independent-specimens.wgswxs.primary-plus.tsv
    ├── independent-specimens.wgswxs.primary.tsv
    ├── intersect_cds_lancet_strelka_mutect_WGS.bed
    ├── intersect_cds_lancet_WXS.bed
    ├── intersect_strelka_mutect_WGS.bed
    ├── md5sum.txt
    ├── pbta-cnv-cnvkit-gistic.zip
    ├── pbta-cnv-cnvkit.seg.gz
    ├── pbta-cnv-consensus-gistic.zip
    ├── pbta-cnv-consensus.seg.gz
    ├── pbta-cnv-controlfreec.tsv.gz
    ├── pbta-fusion-arriba.tsv.gz
    ├── pbta-fusion-putative-oncogenic.tsv
    ├── pbta-fusion-recurrently-fused-genes-byhistology.tsv
    ├── pbta-fusion-recurrently-fused-genes-bysample.tsv
    ├── pbta-fusion-starfusion.tsv.gz
    ├── pbta-gene-counts-rsem-expected_count.polya.rds
    ├── pbta-gene-counts-rsem-expected_count.stranded.rds
    ├── pbta-gene-expression-kallisto.polya.rds
    ├── pbta-gene-expression-kallisto.stranded.rds
    ├── pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
    ├── pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
    ├── pbta-gene-expression-rsem-fpkm.polya.rds
    ├── pbta-gene-expression-rsem-fpkm.stranded.rds
    ├── pbta-gene-expression-rsem-tpm.polya.rds
    ├── pbta-gene-expression-rsem-tpm.stranded.rds
    ├── pbta-histologies.tsv    
    ├── pbta-histologies-base.tsv
    ├── pbta-isoform-counts-rsem-expected_count.polya.rds
    ├── pbta-isoform-counts-rsem-expected_count.stranded.rds
    ├── pbta-isoform-expression-rsem-tpm.polya.rds
    ├── pbta-isoform-expression-rsem-tpm.stranded.rds
    ├── pbta-mend-qc-manifest.tsv
    ├── pbta-mend-qc-results.tar.gz
    ├── pbta-snv-consensus-mutation.maf.tsv.gz
    ├── pbta-snv-consensus-mutation-tmb-all.tsv
    ├── pbta-snv-consensus-mutation-tmb-coding.tsv
    ├── pbta-snv-lancet.vep.maf.gz
    ├── pbta-snv-mutect2.vep.maf.gz
    ├── pbta-snv-strelka2.vep.maf.gz
    ├── pbta-snv-vardict.vep.maf.gz
    ├── pbta-star-log-final.tar.gz
    ├── pbta-star-log-manifest.tsv
    ├── pbta-sv-manta.tsv.gz
    ├── pbta-tcga-manifest.tsv
    ├── pbta-tcga-snv-lancet.vep.maf.gz
    ├── pbta-tcga-snv-mutect2.vep.maf.gz
    ├── pbta-tcga-snv-strelka2.vep.maf.gz
    ├── release-notes.md
    ├── StrexomeLite_hg38_liftover_100bp_padded.bed
    ├── StrexomeLite_Targets_CrossMap_hg38_filtered_chr_prefixed.bed
    ├── WGS.hg38.lancet.300bp_padded.bed
    ├── WGS.hg38.lancet.unpadded.bed
    ├── WGS.hg38.mutect2.vardict.unpadded.bed
    ├── WGS.hg38.strelka2.unpadded.bed
    ├── WGS.hg38.vardict.100bp_padded.bed
    ├── WXS.hg38.100bp_padded.bed
    ├── WXS.hg38.lancet.400bp_padded.bed
    ├── intersected_whole_exome_agilent_designed_120_AND_tcga_6k_genes.Gh38.bed
    ├── intersected_whole_exome_agilent_plus_tcga_6k_AND_tcga_6k_genes.Gh38.bed
    ├── tcga_6k_genes.targetIntervals.Gh38.bed
    ├── tcga_6k_genes.targetIntervals.bed
    ├── tcga-snv-consensus-snv.maf.tsv.gz
    ├── tcga-snv-mutation-tmb-all.tsv
    ├── tcga-snv-mutation-tmb-coding.tsv
    ├── whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.Gh38.bed
    ├── whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.bed
    ├── whole_exome_agilent_designed_120.targetIntervals.Gh38.bed
    ├── whole_exome_agilent_designed_120.targetIntervals.bed
    ├── whole_exome_agilent_plus_tcga_6k.targetIntervals.Gh38.bed
    └── whole_exome_agilent_plus_tcga_6k.targetIntervals.bed
```

## archived release
### release-v17-20200908
- release date: 2020-09-11
- status: available
- changes:
  - pbta-histologies.tsv updates:
    - clinical updates to columns: `aliquot_id`, `parent_aliquot_id`, `reported_gender`, `race`, `ethnicity`, `tumor_descriptor`, `primary_site`, `age_at_diagnosis_days`, `pathology_diagnosis`, `integrated_diagnosis`, `short_histology`, `broad_histology`, `Notes`, `OS_days`, `OS_status`, `age_last_update_days`,`seq_center`, `cancer_predispositions`, `molecular_subtype`
    - Updated logic for `integrated_diagnosis`, `Notes`, `short_histology` and `broad_histology`: If the sample was subtyped, values were updated based on new subtyping, if the `pathology_diagnosis` changed between v16 and v17, all values were set to `NA` and will be updated in the next release, all remaining samples were set to the values from v16.
    - changed column `glioma_brain_region` to `CNS_region`, updated logic and entities, and assigned regions for all tumors
    - added columns `pathology_free_text_diagnosis` and `cohort_participant_id`
    - added sequencing center for `BS_QB84TBA3` and `BS_S2TA8R29` (previously NA)
    - aliquot_id for `BS_MVYA262V` sample_id `7316-3213` was changed from `A05265` to `A05233`
    - updated medulloblastoma subtypes per [#742](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/742)
    - added LGAT subtyping per [#631](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/631)
    - added EWS diagnoses per [#623](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/623)
  - update PBTA TMB file `pbta-snv-consensus-mutation-tmb-coding.tsv` post bug fix, etc. per [#727](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/727) and [#728](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/728)
  - new consensus mutation and TMB files from [TCGA analysis](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers) :
    - tcga-snv-consensus-snv.maf.tsv.gz
    - tcga-snv-mutation-tmb-all.tsv
    - tcga-snv-mutation-tmb-coding.tsv
- folder structure:
```
data
└── release-v17-20200908
    ├── consensus_seg_annotated_cn_autosomes.tsv.gz
    ├── consensus_seg_annotated_cn_x_and_y.tsv.gz
    ├── data-files-description.md
    ├── fusion_summary_embryonal_foi.tsv
    ├── fusion_summary_ependymoma_foi.tsv
    ├── fusion_summary_ewings_foi.tsv
    ├── independent-specimens.wgs.primary-plus.tsv
    ├── independent-specimens.wgs.primary.tsv
    ├── independent-specimens.wgswxs.primary-plus.tsv
    ├── independent-specimens.wgswxs.primary.tsv
    ├── intersect_cds_lancet_strelka_mutect_WGS.bed
    ├── intersect_cds_lancet_WXS.bed
    ├── intersect_strelka_mutect_WGS.bed
    ├── md5sum.txt
    ├── pbta-cnv-cnvkit-gistic.zip
    ├── pbta-cnv-cnvkit.seg.gz
    ├── pbta-cnv-consensus-gistic.zip
    ├── pbta-cnv-consensus.seg.gz
    ├── pbta-cnv-controlfreec.tsv.gz
    ├── pbta-fusion-arriba.tsv.gz
    ├── pbta-fusion-putative-oncogenic.tsv
    ├── pbta-fusion-recurrently-fused-genes-byhistology.tsv
    ├── pbta-fusion-recurrently-fused-genes-bysample.tsv
    ├── pbta-fusion-starfusion.tsv.gz
    ├── pbta-gene-counts-rsem-expected_count.polya.rds
    ├── pbta-gene-counts-rsem-expected_count.stranded.rds
    ├── pbta-gene-expression-kallisto.polya.rds
    ├── pbta-gene-expression-kallisto.stranded.rds
    ├── pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
    ├── pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
    ├── pbta-gene-expression-rsem-fpkm.polya.rds
    ├── pbta-gene-expression-rsem-fpkm.stranded.rds
    ├── pbta-gene-expression-rsem-tpm.polya.rds
    ├── pbta-gene-expression-rsem-tpm.stranded.rds
    ├── pbta-histologies.tsv
    ├── pbta-isoform-counts-rsem-expected_count.polya.rds
    ├── pbta-isoform-counts-rsem-expected_count.stranded.rds
    ├── pbta-isoform-expression-rsem-tpm.polya.rds
    ├── pbta-isoform-expression-rsem-tpm.stranded.rds
    ├── pbta-mend-qc-manifest.tsv
    ├── pbta-mend-qc-results.tar.gz
    ├── pbta-snv-consensus-mutation.maf.tsv.gz
    ├── pbta-snv-consensus-mutation-tmb-all.tsv
    ├── pbta-snv-consensus-mutation-tmb-coding.tsv
    ├── pbta-snv-lancet.vep.maf.gz
    ├── pbta-snv-mutect2.vep.maf.gz
    ├── pbta-snv-strelka2.vep.maf.gz
    ├── pbta-snv-vardict.vep.maf.gz
    ├── pbta-star-log-final.tar.gz
    ├── pbta-star-log-manifest.tsv
    ├── pbta-sv-manta.tsv.gz
    ├── pbta-tcga-manifest.tsv
    ├── pbta-tcga-snv-lancet.vep.maf.gz
    ├── pbta-tcga-snv-mutect2.vep.maf.gz
    ├── pbta-tcga-snv-strelka2.vep.maf.gz
    ├── release-notes.md
    ├── StrexomeLite_hg38_liftover_100bp_padded.bed
    ├── StrexomeLite_Targets_CrossMap_hg38_filtered_chr_prefixed.bed
    ├── WGS.hg38.lancet.300bp_padded.bed
    ├── WGS.hg38.lancet.unpadded.bed
    ├── WGS.hg38.mutect2.vardict.unpadded.bed
    ├── WGS.hg38.strelka2.unpadded.bed
    ├── WGS.hg38.vardict.100bp_padded.bed
    ├── WXS.hg38.100bp_padded.bed
    ├── WXS.hg38.lancet.400bp_padded.bed
    ├── intersected_whole_exome_agilent_designed_120_AND_tcga_6k_genes.Gh38.bed
    ├── intersected_whole_exome_agilent_plus_tcga_6k_AND_tcga_6k_genes.Gh38.bed
    ├── tcga_6k_genes.targetIntervals.Gh38.bed
    ├── tcga_6k_genes.targetIntervals.bed
    ├── tcga-snv-consensus-snv.maf.tsv.gz
    ├── tcga-snv-mutation-tmb-all.tsv
    ├── tcga-snv-mutation-tmb-coding.tsv
    ├── whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.Gh38.bed
    ├── whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.bed
    ├── whole_exome_agilent_designed_120.targetIntervals.Gh38.bed
    ├── whole_exome_agilent_designed_120.targetIntervals.bed
    ├── whole_exome_agilent_plus_tcga_6k.targetIntervals.Gh38.bed
    └── whole_exome_agilent_plus_tcga_6k.targetIntervals.bed
```

##archived release
### release-v16-20200320
- release date: 2020-03-27
- status: available
- changes:
  - Updated TCGA MAFs per #568
  - Update TCGA manifest file to include BED file used per sample
  - Add new TCGA BED files, padded and unpadded
  - Remove old TCGA BED file `gencode.v19.basic.exome.hg38liftover.100bp_padded.bed`
  - new fusions files from #621:
    - `pbta-fusion-putative-oncogenic.tsv`
    - `pbta-fusion-recurrently-fused-genes-byhistology.tsv`
    - `pbta-fusion-recurrently-fused-genes-bysample.tsv`
  - new fusion summarys from `analyses/fusion-summary`
    - `fusion_summary_ewings_foi.tsv`
- folder structure:
```
data
└── release-v16-20200320
    ├── consensus_seg_annotated_cn_autosomes.tsv.gz
    ├── consensus_seg_annotated_cn_x_and_y.tsv.gz
    ├── data-files-description.md
    ├── fusion_summary_embryonal_foi.tsv
    ├── fusion_summary_ependymoma_foi.tsv
    ├── fusion_summary_ewings_foi.tsv
    ├── independent-specimens.wgs.primary-plus.tsv
    ├── independent-specimens.wgs.primary.tsv
    ├── independent-specimens.wgswxs.primary-plus.tsv
    ├── independent-specimens.wgswxs.primary.tsv
    ├── intersect_cds_lancet_strelka_mutect_WGS.bed
    ├── intersect_cds_lancet_WXS.bed
    ├── intersect_strelka_mutect_WGS.bed
    ├── md5sum.txt
    ├── pbta-cnv-cnvkit-gistic.zip
    ├── pbta-cnv-cnvkit.seg.gz
    ├── pbta-cnv-consensus-gistic.zip
    ├── pbta-cnv-consensus.seg.gz
    ├── pbta-cnv-controlfreec.tsv.gz
    ├── pbta-fusion-arriba.tsv.gz
    ├── pbta-fusion-putative-oncogenic.tsv
    ├── pbta-fusion-recurrently-fused-genes-byhistology.tsv
    ├── pbta-fusion-recurrently-fused-genes-bysample.tsv
    ├── pbta-fusion-starfusion.tsv.gz
    ├── pbta-gene-counts-rsem-expected_count.polya.rds
    ├── pbta-gene-counts-rsem-expected_count.stranded.rds
    ├── pbta-gene-expression-kallisto.polya.rds
    ├── pbta-gene-expression-kallisto.stranded.rds
    ├── pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
    ├── pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
    ├── pbta-gene-expression-rsem-fpkm.polya.rds
    ├── pbta-gene-expression-rsem-fpkm.stranded.rds
    ├── pbta-gene-expression-rsem-tpm.polya.rds
    ├── pbta-gene-expression-rsem-tpm.stranded.rds
    ├── pbta-histologies.tsv
    ├── pbta-isoform-counts-rsem-expected_count.polya.rds
    ├── pbta-isoform-counts-rsem-expected_count.stranded.rds
    ├── pbta-isoform-expression-rsem-tpm.polya.rds
    ├── pbta-isoform-expression-rsem-tpm.stranded.rds
    ├── pbta-mend-qc-manifest.tsv
    ├── pbta-mend-qc-results.tar.gz
    ├── pbta-snv-consensus-mutation.maf.tsv.gz
    ├── pbta-snv-consensus-mutation-tmb-all.tsv
    ├── pbta-snv-consensus-mutation-tmb-coding.tsv
    ├── pbta-snv-lancet.vep.maf.gz
    ├── pbta-snv-mutect2.vep.maf.gz
    ├── pbta-snv-strelka2.vep.maf.gz
    ├── pbta-snv-vardict.vep.maf.gz
    ├── pbta-star-log-final.tar.gz
    ├── pbta-star-log-manifest.tsv
    ├── pbta-sv-manta.tsv.gz
    ├── pbta-tcga-manifest.tsv
    ├── pbta-tcga-snv-lancet.vep.maf.gz
    ├── pbta-tcga-snv-mutect2.vep.maf.gz
    ├── pbta-tcga-snv-strelka2.vep.maf.gz
    ├── release-notes.md
    ├── StrexomeLite_hg38_liftover_100bp_padded.bed
    ├── StrexomeLite_Targets_CrossMap_hg38_filtered_chr_prefixed.bed
    ├── WGS.hg38.lancet.300bp_padded.bed
    ├── WGS.hg38.lancet.unpadded.bed
    ├── WGS.hg38.mutect2.vardict.unpadded.bed
    ├── WGS.hg38.strelka2.unpadded.bed
    ├── WGS.hg38.vardict.100bp_padded.bed
    ├── WXS.hg38.100bp_padded.bed
    ├── WXS.hg38.lancet.400bp_padded.bed
    ├── intersected_whole_exome_agilent_designed_120_AND_tcga_6k_genes.Gh38.bed
    ├── intersected_whole_exome_agilent_plus_tcga_6k_AND_tcga_6k_genes.Gh38.bed
    ├── tcga_6k_genes.targetIntervals.Gh38.bed
    ├── tcga_6k_genes.targetIntervals.bed
    ├── whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.Gh38.bed
    ├── whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.bed
    ├── whole_exome_agilent_designed_120.targetIntervals.Gh38.bed
    ├── whole_exome_agilent_designed_120.targetIntervals.bed
    ├── whole_exome_agilent_plus_tcga_6k.targetIntervals.Gh38.bed
    └── whole_exome_agilent_plus_tcga_6k.targetIntervals.bed
```

## archived release
### release-v15-20200228
- release date: 2020-02-28
- status: available
- changes:
  - Update Strelka2, Mutect2, and Lancet TCGA MAF files per [#257](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/257); files updated:
    - pbta-tcga-snv-mutect2.vep.maf.gz
    - pbta-tcga-snv-strelka2.vep.maf.gz
    - pbta-tcga-snv-lancet.vep.maf.gz
  - Add TCGA WXS BED regions file lifted from hg19 to hg38 and padded by 100bp, which was used for variant calling:
    - gencode.v19.basic.exome.hg38liftover.100bp_padded.bed, obtained from [the GDC website](https://gdc.cancer.gov/about-data/publications/mc3-20170); NOTE: this BED file is almost certainly wrong per [#568](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/568).
  - Replace pbta-cnv-consensus-gistic.zip with the file from the [OpenPBTA analysis master repository](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/46cf6ccb119312ccae6122ac94c51710df01f6da/analyses/run-gistic/results) now that GISTIC has been installed [comment](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/543#issuecomment-587514075).
  - Updated `pbta-histologies.tsv` to:
    - Harmonize `Embryonal Tumor` and `Embryonal tumor` [#541](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/541)
    - Change `disease_type_old` to `pathology_diagnosis` and `disease_type_new` to `integrated_diagnosis` per request in [comment](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/543#issuecomment-589248029).
  - Add consensus focal CN files from [analysis](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/46cf6ccb119312ccae6122ac94c51710df01f6da/analyses/focal-cn-file-preparation):
    - consensus_seg_annotated_cn_autosomes.tsv.gz
    - consensus_seg_annotated_cn_x_and_y.tsv.gz
  - Updated fusion files to add additional kinase genes, remove residual polyA stranded samples, andfix order of filtering operations per [comment](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/539#issuecomment-587167746), [comment](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/521#issuecomment-582990375), [#530](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/530), and [#553](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/530), and [PR #567](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/567):
    - pbta-fusion-recurrently-fused-genes-byhistology.tsv
    - pbta-fusion-putative-oncogenic.tsv
    - pbta-fusion-recurrently-fused-genes-bysample.tsv
- folder structure:
```
data
└── release-v15-20200228
    ├── release-notes.md
    ├── data-files-description.md
    ├── StrexomeLite_Targets_CrossMap_hg38_filtered_chr_prefixed.bed
    ├── StrexomeLite_hg38_liftover_100bp_padded.bed
    ├── WGS.hg38.lancet.300bp_padded.bed
    ├── WGS.hg38.lancet.unpadded.bed
    ├── WGS.hg38.mutect2.vardict.unpadded.bed
    ├── WGS.hg38.strelka2.unpadded.bed
    ├── WGS.hg38.vardict.100bp_padded.bed
    ├── WXS.hg38.100bp_padded.bed
    ├── WXS.hg38.lancet.400bp_padded.bed
    ├── md5sum.txt
    ├── pbta-cnv-cnvkit.seg.gz
    ├── pbta-cnv-consensus.seg.gz
    ├── pbta-cnv-controlfreec.tsv.gz
    ├── pbta-cnv-cnvkit-gistic.zip
    ├── pbta-cnv-consensus-gistic.zip
    ├── pbta-fusion-arriba.tsv.gz
    ├── pbta-fusion-starfusion.tsv.gz
    ├── pbta-fusion-putative-oncogenic.tsv
    ├── pbta-gene-counts-rsem-expected_count.polya.rds
    ├── pbta-gene-counts-rsem-expected_count.stranded.rds
    ├── pbta-gene-expression-kallisto.polya.rds
    ├── pbta-gene-expression-kallisto.stranded.rds
    ├── pbta-gene-expression-rsem-fpkm.polya.rds
    ├── pbta-gene-expression-rsem-fpkm.stranded.rds
    ├── pbta-histologies.tsv
    ├── pbta-snv-lancet.vep.maf.gz
    ├── pbta-snv-mutect2.vep.maf.gz
    ├── pbta-snv-strelka2.vep.maf.gz
    ├── pbta-snv-vardict.vep.maf.gz
    ├── pbta-sv-manta.tsv.gz
    ├── independent-specimens.wgs.primary-plus.tsv
    ├── independent-specimens.wgs.primary.tsv
    ├── independent-specimens.wgswxs.primary-plus.tsv
    ├── independent-specimens.wgswxs.primary.tsv
    ├── pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
    ├── pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
    ├── pbta-gene-expression-rsem-tpm.polya.rds
    ├── pbta-gene-expression-rsem-tpm.stranded.rds
    ├── pbta-isoform-expression-rsem-tpm.polya.rds
    ├── pbta-isoform-expression-rsem-tpm.stranded.rds
    ├── pbta-isoform-counts-rsem-expected_count.polya.rds
    ├── pbta-isoform-counts-rsem-expected_count.stranded.rds
    ├── pbta-snv-consensus-mutation.maf.tsv.gz
    ├── pbta-snv-consensus-mutation-tmb-all.tsv
    ├── pbta-snv-consensus-mutation-tmb-coding.tsv
    ├── pbta-fusion-recurrently-fused-genes-byhistology.tsv
    ├── pbta-fusion-recurrently-fused-genes-bysample.tsv
    ├── pbta-tcga-snv-lancet.vep.maf.gz
    ├── pbta-tcga-snv-strelka2.vep.maf.gz
    ├── pbta-tcga-snv-mutect2.vep.maf.gz
    ├── pbta-tcga-manifest.tsv
    ├── pbta-mend-qc-results.tar.gz
    ├── pbta-mend-qc-manifest.tsv
    ├── pbta-star-log-final.tar.gz
    ├── pbta-star-log-manifest.tsv
    ├── intersect_cds_lancet_strelka_mutect_WGS.bed
    ├── intersect_cds_lancet.bed
    ├── intersect_strelka_mutect_WGS.bed
    ├── fusion_summary_embryonal_foi.tsv
    ├── fusion_summary_ependymoma_foi.tsv
    ├── gencode.v19.basic.exome.hg38liftover.100bp_padded.bed
    ├── consensus_seg_annotated_cn_autosomes.tsv.gz
    └── consensus_seg_annotated_cn_x_and_y.tsv.gz
```

### release-v14-20200203
- release date: 2020-02-03
- status: available
- changes:
  - Update kallisto stranded file to remove index column per [#474](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/474). Also removed samples from last re-sequencing polyA+stranded batch that were missed with the v13 release:
    - pbta-gene-expression-kallisto.stranded.rds
  - Update matrices of ependymonal tumor and embryonal tumor fusions of interest by biospecimen from [`analyses/fusion-summary`](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/478) to include all RNA biospecimens in new `pbta-histologies.tsv` file without fusion calls. Files updated:
    - fusion_summary_embryonal_foi.tsv
    - fusion_summary_ependymoma_foi.tsv
  - Update Strelka2, Mutect2, and Lancet TCGA MAF files per [#483](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/483) and [#512](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/512). Files updated:
    - pbta-tcga-snv-mutect2.vep.maf.gz
    - pbta-tcga-snv-strelka2.vep.maf.gz
    - pbta-tcga-snv-lancet.vep.maf.gz
  - Remove `cnv_consensus.tsv` file per [this comment](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/432#issuecomment-579340589).
  - Update copy number files:
    - Update CNVkit seg file to add missing specimen per [#472](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/472):
        - pbta-cnv-cnvkit.seg.gz
    - Update GISTIC results for CNVkit to include missing specimen per [#491](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/491):
        - pbta-cnv-cnvkit-gistic.zip
    - Add consensus SEG file per [#441](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/441):
        - pbta-cnv-consensus.seg.gz
    - Add GISTIC results for consensus SEG per [#453](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/453):
        - pbta-cnv-consensus-gistic.zip
  - Update analysis files and names from `exon` to `cds` per [#440](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/440) and [this comment](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/432#issuecomment-581462059): 
    - intersect_strelka_mutect_WGS.bed
    - intersect_cds_lancet.bed
    - intersect_cds_lancet_strelka_mutect_WGS.bed
    - pbta-snv-consensus-mutation-tmb-coding.tsv
  - Update `pbta-histologies.tsv` to add embryonal `molecular_subtypes` per [#251](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/251) using results [here](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/molecular-subtyping-embryonal/results) and high-grade glioma (HGG) `molecular_subtypes` per [#249](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/249) and [this commit](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/435/commits/7e8e7917c002e11ab97312dc7f64c542fc04892b). 
    - Additional clinical data updated include:
        - `glioma_brain_region` (if sample was not previously classified as glioma)
        - `Notes` (to concatenate old notes for `disease_type_new` and `molecular_subtype` changes and current changes based on OpenPBTA subtyping)
        - `disease_type_new` 
            - HGG: created "Diffuse midline glioma" for any subtype containing `DMG` and rest were "High-grade glioma". (Previously, DMGs were both HGG and DIPG, but following WHO 2016 nomenclature, we stick with DMG instead of DIPG now that we have subtypes).
                - Embryonal: if `molecular_subtype` contains `ETMR`, these became `Embryonal tumor with multilayer rosettes` while the rest became `CNS Embryonal Tumor`        
        - `short_histology`
            - HGG: all became `HGAT`
            - Embryonal: if `molecular_subtype` contains `ETMR`, these became `ETMR` while the rest became `Embryonal Tumor`
        - `broad_histology`
            - HGG: all became `Diffuse astrocytic and oligodendroglial tumor`
            - Embryonal: all became `Embryonal Tumor`

- folder structure:
```
data
└── release-v14-20200203
    ├── release-notes.md
    ├── data-files-description.md
    ├── StrexomeLite_Targets_CrossMap_hg38_filtered_chr_prefixed.bed
    ├── StrexomeLite_hg38_liftover_100bp_padded.bed
    ├── WGS.hg38.lancet.300bp_padded.bed
    ├── WGS.hg38.lancet.unpadded.bed
    ├── WGS.hg38.mutect2.vardict.unpadded.bed
    ├── WGS.hg38.strelka2.unpadded.bed
    ├── WGS.hg38.vardict.100bp_padded.bed
    ├── WXS.hg38.100bp_padded.bed
    ├── WXS.hg38.lancet.400bp_padded.bed
    ├── md5sum.txt
    ├── pbta-cnv-cnvkit.seg.gz
    ├── pbta-cnv-consensus.seg.gz
    ├── pbta-cnv-controlfreec.tsv.gz
    ├── pbta-cnv-cnvkit-gistic.zip
    ├── pbta-cnv-consensus-gistic.zip
    ├── pbta-fusion-arriba.tsv.gz
    ├── pbta-fusion-starfusion.tsv.gz
    ├── pbta-fusion-putative-oncogenic.tsv
    ├── pbta-gene-counts-rsem-expected_count.polya.rds
    ├── pbta-gene-counts-rsem-expected_count.stranded.rds
    ├── pbta-gene-expression-kallisto.polya.rds
    ├── pbta-gene-expression-kallisto.stranded.rds
    ├── pbta-gene-expression-rsem-fpkm.polya.rds
    ├── pbta-gene-expression-rsem-fpkm.stranded.rds
    ├── pbta-histologies.tsv
    ├── pbta-snv-lancet.vep.maf.gz
    ├── pbta-snv-mutect2.vep.maf.gz
    ├── pbta-snv-strelka2.vep.maf.gz
    ├── pbta-snv-vardict.vep.maf.gz
    ├── pbta-sv-manta.tsv.gz
    ├── independent-specimens.wgs.primary-plus.tsv
    ├── independent-specimens.wgs.primary.tsv
    ├── independent-specimens.wgswxs.primary-plus.tsv
    ├── independent-specimens.wgswxs.primary.tsv
    ├── pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
    ├── pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
    ├── pbta-gene-expression-rsem-tpm.polya.rds
    ├── pbta-gene-expression-rsem-tpm.stranded.rds
    ├── pbta-isoform-expression-rsem-tpm.polya.rds
    ├── pbta-isoform-expression-rsem-tpm.stranded.rds
    ├── pbta-isoform-counts-rsem-expected_count.polya.rds
    ├── pbta-isoform-counts-rsem-expected_count.stranded.rds
    ├── pbta-snv-consensus-mutation.maf.tsv.gz
    ├── pbta-snv-consensus-mutation-tmb-all.tsv
    ├── pbta-snv-consensus-mutation-tmb-coding.tsv
    ├── pbta-fusion-recurrently-fused-genes-byhistology.tsv
    ├── pbta-fusion-recurrently-fused-genes-bysample.tsv
    ├── pbta-tcga-snv-lancet.vep.maf.gz
    ├── pbta-tcga-snv-strelka2.vep.maf.gz
    ├── pbta-tcga-snv-mutect2.vep.maf.gz
    ├── pbta-tcga-manifest.tsv
    ├── pbta-mend-qc-results.tar.gz
    ├── pbta-mend-qc-manifest.tsv
    ├── pbta-star-log-final.tar.gz
    ├── pbta-star-log-manifest.tsv
    ├── intersect_cds_lancet_strelka_mutect_WGS.bed
    ├── intersect_cds_lancet.bed
    ├── intersect_strelka_mutect_WGS.bed
    ├── fusion_summary_embryonal_foi.tsv
    └── fusion_summary_ependymoma_foi.tsv
```

## archived releases
### release-v13-20200116
- release date: 2020-01-16
- status: available
- changes:
  - Remove stranded RNA-Seq for 23 PNOC samples and 21 CBTTC samples previously sequenced using a polyA library prep and simultaneously address issues [#369](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/369), [#370](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/370), and [#371](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/369). Files updated:
    - pbta-fusion-arriba.tsv.gz
    - pbta-fusion-starfusion.tsv.gz
    - pbta-gene-expression-rsem-tpm.stranded.rds
    - pbta-gene-expression-rsem-fpkm.stranded.rds
    - pbta-isoform-expression-rsem-tpm.stranded.rds
    - pbta-isoform-counts-rsem-expected_count.stranded.rds
    - pbta-gene-counts-rsem-expected_count.stranded.rds
    - pbta-gene-expression-kallisto.stranded.rds
    - pbta-fusion-recurrently-fused-genes-byhistology.tsv
    - pbta-fusion-recurrently-fused-genes-bysample.tsv
    - pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
    - pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
  - Add matrices of ependymonal tumor and embryonal tumor fusions of interest by biospecimen from [`analyses/fusion-summary`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/fusion-summary), subsetted for RNA biospecimens in new `pbta-histologies.tsv` file. Files added:
    - fusion_summary_embryonal_foi.tsv
    - fusion_summary_ependymoma_foi.tsv
  - Add MendQC and STAR output files and associated manifests, per [#341](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/341). Files added:
    - pbta-mend-qc-results.tar.gz
    - pbta-mend-qc-manifest.tsv
    - pbta-star-log-final.tar.gz
    - pbta-star-log-manifest.tsv
  - Add TCGA MAF and clinical files per [#257](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/257). Files added:
    - pbta-tcga-snv-lancet.vep.maf.gz
    - pbta-tcga-snv-mutect2.vep.maf.gz
    - pbta-tcga-snv-strelka2.vep.maf.gz
    - pbta-tcga-manifest.tsv
  - Add `cnv_consensus.tsv` file from [#128](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/128).
  - Add analysis files from [#351](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/351#issuecomment-571295623)
    - intersect_strelka_mutect_WGS.bed
    - intersect_exon_lancet_strelka_mutect_WGS.bed
    - intersect_exon_WXS.bed
  - Add `WXS.hg38.lancet.400bp_padded.bed` file. 
  - Update `pbta-histologies.tsv` to remove RNA-Seq samples listed above, propagate medulloblastoma `molecular_subtypes` per [#379](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/379), harmonize "Diagnosis" to "Initial CNS Tumor", fix PNOC003 `seq_center` for RNA-Seq samples to "TGEN", harmonize `ethnicity` to match CBTTC data in [KidsFirst](https://kidsfirstdrc.org/):


 | Old ethnicity | New ethnicity |
 |---------------|--------------------------|
 | Non-hispanic | Not Hispanic or Latino |
 | Unknown | Unavailable/Not Reported |
 | Not Reported | Unavailable/Not Reported |


- folder structure:
```
data
└── release-v13-20200116
    ├── release-notes.md
    ├── data-files-description.md
    ├── StrexomeLite_Targets_CrossMap_hg38_filtered_chr_prefixed.bed
    ├── StrexomeLite_hg38_liftover_100bp_padded.bed
    ├── WGS.hg38.lancet.300bp_padded.bed
    ├── WGS.hg38.lancet.unpadded.bed
    ├── WGS.hg38.mutect2.vardict.unpadded.bed
    ├── WGS.hg38.strelka2.unpadded.bed
    ├── WGS.hg38.vardict.100bp_padded.bed
    ├── WXS.hg38.100bp_padded.bed
    ├── WXS.hg38.lancet.400bp_padded.bed
    ├── md5sum.txt
    ├── pbta-cnv-cnvkit.seg.gz
    ├── pbta-cnv-controlfreec.tsv.gz
    ├── pbta-cnv-cnvkit-gistic.zip
    ├── pbta-fusion-arriba.tsv.gz
    ├── pbta-fusion-starfusion.tsv.gz
    ├── pbta-fusion-putative-oncogenic.tsv
    ├── pbta-gene-counts-rsem-expected_count.polya.rds
    ├── pbta-gene-counts-rsem-expected_count.stranded.rds
    ├── pbta-gene-expression-kallisto.polya.rds
    ├── pbta-gene-expression-kallisto.stranded.rds
    ├── pbta-gene-expression-rsem-fpkm.polya.rds
    ├── pbta-gene-expression-rsem-fpkm.stranded.rds
    ├── pbta-histologies.tsv
    ├── pbta-snv-lancet.vep.maf.gz
    ├── pbta-snv-mutect2.vep.maf.gz
    ├── pbta-snv-strelka2.vep.maf.gz
    ├── pbta-snv-vardict.vep.maf.gz
    ├── pbta-sv-manta.tsv.gz
    ├── independent-specimens.wgs.primary-plus.tsv
    ├── independent-specimens.wgs.primary.tsv
    ├── independent-specimens.wgswxs.primary-plus.tsv
    ├── independent-specimens.wgswxs.primary.tsv
    ├── pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
    ├── pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
    ├── pbta-gene-expression-rsem-tpm.polya.rds
    ├── pbta-gene-expression-rsem-tpm.stranded.rds
    ├── pbta-isoform-expression-rsem-tpm.polya.rds
    ├── pbta-isoform-expression-rsem-tpm.stranded.rds
    ├── pbta-isoform-counts-rsem-expected_count.polya.rds
    ├── pbta-isoform-counts-rsem-expected_count.stranded.rds
    ├── pbta-snv-consensus-mutation.maf.tsv.gz
    ├── pbta-snv-consensus-mutation-tmb-all.tsv
    ├── pbta-snv-consensus-mutation-tmb-coding.tsv
    ├── pbta-fusion-recurrently-fused-genes-byhistology.tsv
    ├── pbta-fusion-recurrently-fused-genes-bysample.tsv
    ├── pbta-tcga-snv-lancet.vep.maf.gz
    ├── pbta-tcga-snv-strelka2.vep.maf.gz
    ├── pbta-tcga-snv-mutect2.vep.maf.gz
    ├── pbta-tcga-manifest.tsv
    ├── pbta-mend-qc-results.tar.gz
    ├── pbta-mend-qc-manifest.tsv
    ├── pbta-star-log-final.tar.gz
    ├── pbta-star-log-manifest.tsv
    ├── cnv_consensus.tsv
    ├── intersect_exon_lancet_strelka_mutect_WGS.bed
    ├── intersect_exon_WXS.bed
    ├── intersect_strelka_mutect_WGS.bed
    ├── fusion_summary_embryonal_foi.tsv
    └── fusion_summary_ependymoma_foi.tsv
```

### release-v12-20191217
- release date: 2019-12-17
- status: available
- changes:
  - Add `data-file-descriptions.md` with data release to better track file types, origins, and workflows per [#334](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/334) and [#336](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/336)
  - Add stranded RNA-Seq for 23 PNOC samples and 21 CBTTC samples previously sequenced using a polyA library prep. Files updated:
    - pbta-fusion-arriba.tsv.gz
    - pbta-fusion-starfusion.tsv.gz
    - pbta-gene-expression-rsem-tpm.stranded.rds
    - pbta-gene-expression-rsem-fpkm.stranded.rds
    - pbta-isoform-expression-rsem-tpm.stranded.rds
    - pbta-isoform-counts-rsem-expected_count.stranded.rds
    - pbta-gene-counts-rsem-expected_count.stranded.rds
    - pbta-gene-expression-kallisto.stranded.rds
    - pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
  - Add recurrently-fused genes by histology and matrix of recurrently-fused genes by biospecimen from [fusion filtering and prioritization analysis](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/fusion_filtering)
  - Update consensus TMB files and MAF [#333](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/333)
  - Add RNA-Seq [collapsed matrices](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/287) - wrong files (tables of transcripts removed) were included with [V10](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/273)
  - Rename `WGS.hg38.mutect2.unpadded.bed` to `WGS.hg38.mutect2.vardict.unpadded.bed` to better reflect usage
  - Update `pbta-histologies.tsv` to add new RNA-Seq samples listed above, [#222 harmonize disease separators](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/222), and reran [medulloblastoma classifier](https://github.com/d3b-center/medullo-classifier-package) using V12 RSEM fpkm collapsed files
    - BS_2Z1MKS84, BS_5VQP0E6K re-classified from Group4 to WNT and BS_3BDAG9YN, BS_8T7DZV2F, and BS_JTMXAMB7 re-classified from Group3 to WNT
  - Add CNVkit GISTIC results focal CN analyses, eg: [#244](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/244) and [#8](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/8)
- folder structure:
```
data
└── release-v12-20191217
    ├── CHANGELOG.md
    ├── data-file-descriptions.md
    ├── StrexomeLite_Targets_CrossMap_hg38_filtered_chr_prefixed.bed
    ├── StrexomeLite_hg38_liftover_100bp_padded.bed
    ├── WGS.hg38.lancet.300bp_padded.bed
    ├── WGS.hg38.lancet.unpadded.bed
    ├── WGS.hg38.mutect2.vardict.unpadded.bed
    ├── WGS.hg38.strelka2.unpadded.bed
    ├── WGS.hg38.vardict.100bp_padded.bed
    ├── WXS.hg38.100bp_padded.bed
    ├── md5sum.txt
    ├── pbta-cnv-cnvkit.seg.gz
    ├── pbta-cnv-controlfreec.tsv.gz
    ├── pbta-fusion-arriba.tsv.gz
    ├── pbta-fusion-starfusion.tsv.gz
    ├── pbta-fusion-putative-oncogenic.tsv
    ├── pbta-gene-counts-rsem-expected_count.polya.rds
    ├── pbta-gene-counts-rsem-expected_count.stranded.rds
    ├── pbta-gene-expression-kallisto.polya.rds
    ├── pbta-gene-expression-kallisto.stranded.rds
    ├── pbta-gene-expression-rsem-fpkm.polya.rds
    ├── pbta-gene-expression-rsem-fpkm.stranded.rds
    ├── pbta-histologies.tsv
    ├── pbta-isoform-counts-rsem-expected_count.polya.rds
    ├── pbta-isoform-counts-rsem-expected_count.stranded.rds
    ├── pbta-snv-lancet.vep.maf.gz
    ├── pbta-snv-mutect2.vep.maf.gz
    ├── pbta-snv-strelka2.vep.maf.gz
    ├── pbta-snv-vardict.vep.maf.gz
    ├── pbta-sv-manta.tsv.gz
    ├── independent-specimens.wgs.primary-plus.tsv
    ├── independent-specimens.wgs.primary.tsv
    ├── independent-specimens.wgswxs.primary-plus.tsv
    ├── independent-specimens.wgswxs.primary.tsv
    ├── pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
    ├── pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
    ├── pbta-gene-expression-rsem-tpm.polya.rds
    ├── pbta-gene-expression-rsem-tpm.stranded.rds
    ├── pbta-isoform-expression-rsem-tpm.polya.rds
    ├── pbta-isoform-expression-rsem-tpm.stranded.rds
    ├── pbta-snv-consensus-mutation.maf.tsv.gz
    ├── pbta-snv-consensus-mutation-tmb-all.tsv
    ├── pbta-snv-consensus-mutation-tmb-coding.tsv
    ├── pbta-fusion-recurrently-fused-genes-byhistology.tsv
    └── pbta-fusion-recurrently-fused-genes-bysample.tsv
```

### release-v11-20191126
- release date: 2019-11-26
- status: available
- changes:
  - Add putative oncogenic fusions from [fusion filtering and prioritization analysis](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/fusion_filtering)
  - Update consensus SNV/indels to contain [consensus MNVs](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/275)
  - Add RNA-Seq [collapsed matrices](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/287) - wrong files (tables of transcripts removed) were included with [V10](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/273)
- folder structure:
```
data
└── release-v11-20191126
    ├── CHANGELOG.md
    ├── StrexomeLite_Targets_CrossMap_hg38_filtered_chr_prefixed.bed
    ├── StrexomeLite_hg38_liftover_100bp_padded.bed
    ├── WGS.hg38.lancet.300bp_padded.bed
    ├── WGS.hg38.lancet.unpadded.bed
    ├── WGS.hg38.mutect2.unpadded.bed
    ├── WGS.hg38.strelka2.unpadded.bed
    ├── WGS.hg38.vardict.100bp_padded.bed
    ├── WXS.hg38.100bp_padded.bed
    ├── md5sum.txt
    ├── pbta-cnv-cnvkit.seg.gz
    ├── pbta-cnv-controlfreec.tsv.gz
    ├── pbta-fusion-arriba.tsv.gz
    ├── pbta-fusion-starfusion.tsv.gz
    ├── pbta-fusion-putative-oncogenic.tsv
    ├── pbta-gene-counts-rsem-expected_count.polya.rds
    ├── pbta-gene-counts-rsem-expected_count.stranded.rds
    ├── pbta-gene-expression-kallisto.polya.rds
    ├── pbta-gene-expression-kallisto.stranded.rds
    ├── pbta-gene-expression-rsem-fpkm.polya.rds
    ├── pbta-gene-expression-rsem-fpkm.stranded.rds
    ├── pbta-histologies.tsv
    ├── pbta-isoform-counts-rsem-expected_count.polya.rds
    ├── pbta-isoform-counts-rsem-expected_count.stranded.rds
    ├── pbta-snv-lancet.vep.maf.gz
    ├── pbta-snv-mutect2.vep.maf.gz
    ├── pbta-snv-strelka2.vep.maf.gz
    ├── pbta-snv-vardict.vep.maf.gz
    ├── pbta-sv-manta.tsv.gz
    ├── independent-specimens.wgs.primary-plus.tsv
    ├── independent-specimens.wgs.primary.tsv
    ├── independent-specimens.wgswxs.primary-plus.tsv
    ├── independent-specimens.wgswxs.primary.tsv
    ├── pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
    ├── pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
    ├── pbta-gene-expression-rsem-tpm.polya.rds
    ├── pbta-gene-expression-rsem-tpm.stranded.rds
    ├── pbta-isoform-expression-rsem-tpm.polya.rds
    ├── pbta-isoform-expression-rsem-tpm.stranded.rds
    ├── pbta-snv-consensus-mutation.maf.tsv.gz
    └── pbta-snv-consensus-mutation-tmb.tsv
```

### release-v10-20191115
- release date: 2019-11-15
- status: available
- changes:
  - Add RNA-Seq GTF and fasta files per ticket [here](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/241)
  - Add [RSEM gene TPM and isoform matrices](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/229)
  - Add [SNV consensus files](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/207)
  - Add [new MAFs for lancet, vardict, mutect2](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/255)
    -reran VCF2MAF to harmonize columns
  - Add [Collapsed RNA matrices](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/248)
- folder structure:
```
data
└── release-v10-20191115
    ├── CHANGELOG.md
    ├── StrexomeLite_Targets_CrossMap_hg38_filtered_chr_prefixed.bed
    ├── StrexomeLite_hg38_liftover_100bp_padded.bed
    ├── WGS.hg38.lancet.300bp_padded.bed
    ├── WGS.hg38.lancet.unpadded.bed
    ├── WGS.hg38.mutect2.unpadded.bed
    ├── WGS.hg38.strelka2.unpadded.bed
    ├── WGS.hg38.vardict.100bp_padded.bed
    ├── WXS.hg38.100bp_padded.bed
    ├── md5sum.txt
    ├── pbta-cnv-cnvkit.seg.gz
    ├── pbta-cnv-controlfreec.tsv.gz
    ├── pbta-fusion-arriba.tsv.gz
    ├── pbta-fusion-starfusion.tsv.gz
    ├── pbta-gene-counts-rsem-expected_count.polya.rds
    ├── pbta-gene-counts-rsem-expected_count.stranded.rds
    ├── pbta-gene-expression-kallisto.polya.rds
    ├── pbta-gene-expression-kallisto.stranded.rds
    ├── pbta-gene-expression-rsem-fpkm.polya.rds
    ├── pbta-gene-expression-rsem-fpkm.stranded.rds
    ├── pbta-histologies.tsv
    ├── pbta-isoform-counts-rsem-expected_count.polya.rds
    ├── pbta-isoform-counts-rsem-expected_count.stranded.rds
    ├── pbta-snv-lancet.vep.maf.gz
    ├── pbta-snv-mutect2.vep.maf.gz
    ├── pbta-snv-strelka2.vep.maf.gz
    ├── pbta-snv-vardict.vep.maf.gz
    ├── pbta-sv-manta.tsv.gz
    ├── independent-specimens.wgs.primary-plus.tsv
    ├── independent-specimens.wgs.primary.tsv
    ├── independent-specimens.wgswxs.primary-plus.tsv
    ├── independent-specimens.wgswxs.primary.tsv
    ├── pbta-gene-expression-rsem-fpkm-collapsed_table.polya.rds
    ├── pbta-gene-expression-rsem-fpkm-collapsed_table.stranded.rds
    ├── pbta-gene-expression-rsem-tpm.polya.rds
    ├── pbta-gene-expression-rsem-tpm.stranded.rds
    ├── pbta-isoform-expression-rsem-tpm.polya.rds
    ├── pbta-isoform-expression-rsem-tpm.stranded.rds
    └── pbta-snv-consensus_11122019.zip

```

### release-v9-20191105
- release date: 2019-11-05
- status: available
- changes:
  - Updated RNA-Seq FPKM merge files
    - added geneIDs that were missed in previous release
- folder structure:
```
data
└── release-v9-20191105
    ├── CHANGELOG.md
    ├── StrexomeLite_Targets_CrossMap_hg38_filtered_chr_prefixed.bed
    ├── StrexomeLite_hg38_liftover_100bp_padded.bed
    ├── WGS.hg38.lancet.300bp_padded.bed
    ├── WGS.hg38.lancet.unpadded.bed
    ├── WGS.hg38.mutect2.unpadded.bed
    ├── WGS.hg38.strelka2.unpadded.bed
    ├── WGS.hg38.vardict.100bp_padded.bed
    ├── WXS.hg38.100bp_padded.bed
    ├── md5sum.txt
    ├── pbta-cnv-cnvkit.seg.gz
    ├── pbta-cnv-controlfreec.tsv.gz
    ├── pbta-fusion-arriba.tsv.gz
    ├── pbta-fusion-starfusion.tsv.gz
    ├── pbta-gene-counts-rsem-expected_count.polya.rds
    ├── pbta-gene-counts-rsem-expected_count.stranded.rds
    ├── pbta-gene-expression-kallisto.polya.rds
    ├── pbta-gene-expression-kallisto.stranded.rds
    ├── pbta-gene-expression-rsem-fpkm.polya.rds
    ├── pbta-gene-expression-rsem-fpkm.stranded.rds
    ├── pbta-histologies.tsv
    ├── pbta-isoform-counts-rsem-expected_count.polya.rds
    ├── pbta-isoform-counts-rsem-expected_count.stranded.rds
    ├── pbta-snv-lancet.vep.maf.gz
    ├── pbta-snv-mutect2.vep.maf.gz
    ├── pbta-snv-strelka2.vep.maf.gz
    ├── pbta-snv-vardict.vep.maf.gz
    ├── pbta-sv-manta.tsv.gz
    ├── independent-specimens.wgs.primary-plus.tsv
    ├── independent-specimens.wgs.primary.tsv
    ├── independent-specimens.wgswxs.primary-plus.tsv
    └── independent-specimens.wgswxs.primary.tsv
```

### release-v8-20191104
- release date: 2019-11-04
- status: available
- changes:
  - Updated clinical file
    - fixed error in `primary_site`: [#214](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/214)
  - Updated ControlFreeC TSV file
    - added `tumor_ploidy` and updated `genotype` to `segment_genotype`: [PR comment](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/216#discussion_r341868007)
  - Updated RNA-Seq FPKM merge files
    - fixed ID merge error: [#221](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/221) 
- folder structure: same to release-v9-20191105

### release-v7-20191031
- release date: 2019-10-31 :jack_o_lantern:
- status: available
- changes:
  - update molecular_subtype for the clinical file
  - Add ControlFreeC CNV merged TSV file. [data format](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/format/controlfreec-tsv.md)
  - Add ControlFreeC ploidy to clinical file.
  - Add sample lists for analyses requiring unique patient specimens. [issues/155](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/155)
    - independent-specimens.wgs.primary-plus.tsv
    - independent-specimens.wgswxs.primary-plus.tsv
    - independent-specimens.wgs.primary.tsv
    - independent-specimens.wgswxs.primary.tsv

- folder structure: same to release-v9-20191105

### release-v6-20191030
- release date: 2019-10-30
- status: available
- changes:
  - Clinical file updates:
    - Missing `aliquot_id` and `sample_id` added
    - Updated `broad_composition` to `cell line` for WGS samples denoted as Cell line
    - Removed duplicate `BS_4M0ZMCDC` with wrong age at diagnosis
    - Add `cohort` column for CBTTC or PNOC003 samples
    - Tumor specimens missing `composition` were changed to `Solid Tumor`
    - Blood specimens missing `primary_site` were changed to `Peripheral Whole Blood`
    - Updated `age_at_diagnosis` to earliest age reported (same age used in OS calculations)
    - Updated `OS_days` and `OS_status` based on updated clinical data
    - Added `cancer_predispositions` information
    - Added `seq_center` (could not add seq_instrument at this time due to multiple entries for BS_IDs)
    - Harmonized `Diagnosis` and `Initial CNS Tumor` for `tumor_descriptor` field
    - Changed `Relapse` sample to `Progressive` (DIPG sample truly progressive, not relapse)
    - Add tumor purity derived from Theta2 (`normal_fraction` and `tumor_fraction`)
    - Add `glioma_brain_region` for low- and high-grade gliomas
    ​
  - SV:
    - Removed LUMPY data, as additional benchmarking to remove normal SVs needs to be done. We may not include this in a future release.
    ​
  - SNV:
    - Re-ran BS_7KR13R3P using targeted panel bed files; removed WXS calls from MAFs
    - Added WXS calls to all MAFs
    - Added targeted panel bed and padded bed files
    ​
  - CNV:
    - Re-ran ControlFreeC and CNVkit with optional BAF inputs; Added Theta2 purity correction to CNVkit
    - Added copy number to CNVkit and removed ControlFreeC seg file
- folder structure:
```
data
└── release-v6-20191030
    ├── CHANGELOG.md
    ├── StrexomeLite_Targets_CrossMap_hg38_filtered_chr_prefixed.bed
    ├── StrexomeLite_hg38_liftover_100bp_padded.bed
    ├── WGS.hg38.lancet.300bp_padded.bed
    ├── WGS.hg38.lancet.unpadded.bed
    ├── WGS.hg38.mutect2.unpadded.bed
    ├── WGS.hg38.strelka2.unpadded.bed
    ├── WGS.hg38.vardict.100bp_padded.bed
    ├── WXS.hg38.100bp_padded.bed
    ├── md5sum.txt
    ├── pbta-cnv-cnvkit.seg.gz
    ├── pbta-fusion-arriba.tsv.gz
    ├── pbta-fusion-starfusion.tsv.gz
    ├── pbta-gene-counts-rsem-expected_count.polya.rds
    ├── pbta-gene-counts-rsem-expected_count.stranded.rds
    ├── pbta-gene-expression-kallisto.polya.rds
    ├── pbta-gene-expression-kallisto.stranded.rds
    ├── pbta-gene-expression-rsem-fpkm.polya.rds
    ├── pbta-gene-expression-rsem-fpkm.stranded.rds
    ├── pbta-histologies.tsv
    ├── pbta-isoform-counts-rsem-expected_count.polya.rds
    ├── pbta-isoform-counts-rsem-expected_count.stranded.rds
    ├── pbta-snv-lancet.vep.maf.gz
    ├── pbta-snv-mutect2.vep.maf.gz
    ├── pbta-snv-strelka2.vep.maf.gz
    ├── pbta-snv-vardict.vep.maf.gz
    └── pbta-sv-manta.tsv.gz
```

### release-v5-20190924
- release date: 2019-09-24
- status: available
- changes:
  - [Separated RNA-Seq files](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/121):
    - Created separate RDS files for stranded and polyA RNA-Seq samples
  - [new RNA-Seq counts files](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/14):
    - Added RSEM count matrices for genes and transcripts
  - [new ARRIBA file](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/92#discussion_r324873300):
    - Add `annots` column header which was removed during FusionAnnotator run
  - [new SNV files](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/114):
    - Add Lancet VEP-annotated MAF
    - Add VarDict VEP-annotated MAF
  - [new BED interval files](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/3):
    - Add WXS BED (same file used for each variant caller)
    - Add WGS BED files for each variant caller
    - Methods described [here](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md).

- folder structure:
```
data
└── release-v5-20190924
    ├── CHANGELOG.md
    ├── md5sum.txt
    ├── WGS.hg38.lancet.300bp_padded.bed
    ├── WGS.hg38.lancet.unpadded.bed
    ├── WGS.hg38.mutect2.unpadded.bed
    ├── WGS.hg38.strelka2.unpadded.bed
    ├── WGS.hg38.vardict.100bp_padded.bed
    ├── WXS.hg38.100bp_padded.bed
    ├── pbta-cnv-cnvkit.seg.gz
    ├── pbta-cnv-controlfreec.seg.gz
    ├── pbta-fusion-arriba.tsv.gz
    ├── pbta-fusion-starfusion.tsv.gz
    ├── pbta-histologies.tsv
    ├── pbta-snv-mutect2.vep.maf.gz
    ├── pbta-snv-strelka2.vep.maf.gz
    ├── pbta-sv-lumpy.tsv.gz
    ├── pbta-sv-manta.tsv.gz
    ├── pbta-gene-expression-kallisto.polya.rds
    ├── pbta-gene-expression-kallisto.stranded.rds
    ├── pbta-gene-expression-rsem-fpkm.polya.rds
    ├── pbta-gene-expression-rsem-fpkm.stranded.rds
    ├── pbta-gene-counts-rsem-expected_count.polya.rds
    ├── pbta-gene-counts-rsem-expected_count.stranded.rds
    ├── pbta-isoform-counts-rsem-expected_count.polya.rds
    ├── pbta-isoform-counts-rsem-expected_count.stranded.rds
    ├── pbta-snv-lancet.vep.maf.gz
    └── pbta-snv-vardict.vep.maf.gz
```

### release-v4-20190909
- release date: 2019-09-10
- status: available
- changes:
  - [Clinical V4](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/95):
    - De-duplicated `BS_6DCSD5Y6` in the clinical file by removing row with wrong age
    - Added units to `age_at_diagnosis` column (age_at_diagnosis_days)
  - [new RSEM file](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/93):
    - Remove QC-failed samples from V3 RSEM. 
    - Add `BS_XM1AHBDJ` which should be included, but was missing from the V3 RSEM
  - [new LUMPY file](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/97):
    - Add data for 628 more subjects which was missing from V3 LUMPY
- folder structure:
```
data
└── release-v4-20190909
    ├── CHANGELOG.md
    ├── md5sum.txt
    ├── pbta-cnv-cnvkit.seg.gz
    ├── pbta-cnv-controlfreec.seg.gz
    ├── pbta-fusion-arriba.tsv.gz
    ├── pbta-fusion-starfusion.tsv.gz
    ├── pbta-gene-expression-kallisto.rds
    ├── pbta-gene-expression-rsem.fpkm.rds
    ├── pbta-histologies.tsv
    ├── pbta-snv-mutect2.vep.maf.gz
    ├── pbta-snv-strelka2.vep.maf.gz
    ├── pbta-sv-lumpy.tsv.gz
    ├── pbta-sv-manta.tsv.gz
    └── README.md
```

## archived release
### release-v3-20190829
- release date: 2019-08-29
- status: available
- changes:
  - Arriba fusions run with https://github.com/FusionAnnotator/CTAT_HumanFusionLib/wiki#red-herrings-fusion-pairs-that-may-not-be-relevant-to-cancer-and-potential-false-positives to match STAR-fusion results
  - Clinical data sheet, added information for:
    - germline_sex_estimate
    - overall survival (days)
    - overall survival status (living/deceased)
    - RNA_library type
    - edited previous sex variable to be reported_gender
- folder structure: same to `release-v4-20190909`

### release-v2-20190809
- release date: 2019-08-09
- status: available
- changes:
  - removed QC failed RNA-seqs
    - alignment rate < 60%, n=
    - wrong strandedness samples, n=
  - removed tumor/normal genotype mis-match samples
  - removed consent failed samples
  - clincial infomation updated
    - added and harmonized sample ID and diagnosis for PNOC003 samples
    - added medullo subtype information
  - CNV
    - added seg files for CNV data
    - added CNVkit results  
  - SV
    - added LUMPY results
    - added annotated SV files
- folder structure:
```
data
└── release-v2-20190809
    ├── CHANGELOG.md
    ├── md5sum.txt
    ├── pbta-cnv-cnvkit.seg.gz
    ├── pbta-cnv-controlfreec.seg.gz
    ├── pbta-fusion-arriba.tsv.gz
    ├── pbta-fusion-starfusion.tsv.gz
    ├── pbta-gene-expression-kallisto.rds
    ├── pbta-gene-expression-rsem.fpkm.rds
    ├── pbta-histologies.tsv
    ├── pbta-snv-mutect2.vep.maf.gz
    ├── pbta-snv-strelka2.vep.maf.gz
    ├── pbta-sv-lumpy.tsv.gz
    ├── pbta-sv-manta.tsv.gz
    └── README.md
```

### release-v1-20190730
- release date: 2019-07-30
- status: not-available as it contains consent/qc failed samples
- changes: initial push
```
data
└── release-v1-20190730
    ├── arriba.fusions.tsv.gz
    ├── clinical.tsv
    ├── controlfreec.cnv.tsv.gz
    ├── kallisto.abundance.tsv.gz
    ├── kallisto.genes.list
    ├── manta-sv.tsv.gz
    ├── mutect2.maf.gz
    ├── rsem.genes.list
    ├── rsem.genes.tsv.gz
    ├── rsem.isoforms.tsv.gz
    ├── star-fusion.fusions.tsv.gz
    ├── strelka2.maf.gz
    └── tumor-normal-pair.tsv
```
