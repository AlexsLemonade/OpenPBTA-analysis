# release notes
## current release
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
- folder structure:
```
data
└── release-v8-20191104
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

## archived release
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

- folder structure:
```
data
└── release-v7-20191031
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




