# release notes
## current release
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
