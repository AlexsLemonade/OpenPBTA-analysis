# OpenPBTA Open Access Data
This dataset includes all the PBTA somatic mutational and gene expression results in combined tsv or matrix format. We will be releasing this dataset on both [CAVATICA](https://cavatica.sbgenomics.com) and AWS cloud platform moving forward.

## Data format
- CNV results: [SEG format](https://software.broadinstitute.org/software/igv/SEG)
- SNV results: [Annotated MAF format](format/vep-maf.md)
- SV results: [Annotated Manta TSV](format/manta-tsv-header.md)
- Gene Expressions: gene and sample matrix
- Gene Fusions:
    - [Arriba TSV](format/arriba-tsv-header.md)
    - [STARFusion TSV](format/starfusion-tsv-header.md)

## Data Access via CAVATICA
For any user registered on CAVATICA, the latest release of OpenPBTA data can be accessed from the CAVATICA public projects below:
- [Pediatric Brain Tumor Atlas Open Access Data - CBTTC](https://cavatica.sbgenomics.com/u/cavatica/pbta-cbttc/)
- [Pediatric Brain Tumor Atlas Open Access Data - PNOC003](https://cavatica.sbgenomics.com/u/cavatica/pbta-pnoc003/)

## Accessing Data via AWS S3
For any other users, OpenPBTA Open Access Data is also organized by a directory structure for each data release under AWS S3 bucket `s3://kf-openaccess-us-east-1-prd-pbta/data/`. Data folders are named following the naming convention of `release-{version}-{date}`.

Example of the data directory structure:
```
data
└── release-v2-20190809
    ├── release-notes.md
    ├── md5sum.txt
    ├── pbta-cnv-cnvkit.seg.gz
    ├── pbta-cnv-controlfreec.seg.gz
    ├── pbta-fusion-arriba.tsv.gz
    ├── pbta-fusion-starfusion.tsv.gz
    ├── pbta-gene-expression-kallisto.rds
    ├── pbta-histologies.tsv
    ├── pbta-snv-mutect2.vep.maf.gz
    ├── pbta-snv-strelka2.vep.maf.gz
    ├── pbta-sv-manta.tsv.gz
    └── readme.md
```

For the current availabe data, please refer to [release-notes.md](./release-notes.md)

To download data from S3, [aws-cli](https://github.com/aws/aws-cli) is recommened. But other tools like `wget` or `curl` should also work.

Example of download one entire release
```
## aws-cli
aws s3 sync s3://kf-openaccess-us-east-1-prd-pbta/data/release-v2-20190809/ /path/to/local/folder/

## wget
wget https://s3.amazonaws.com/kf-openaccess-us-east-1-prd-pbta/data/release-v2-20190809/md5sum.txt
wget https://s3.amazonaws.com/kf-openaccess-us-east-1-prd-pbta/data/release-v2-20190809/pbta-cnv-cnvkit.seg.gz
wget https://s3.amazonaws.com/kf-openaccess-us-east-1-prd-pbta/data/release-v2-20190809/pbta-cnv-controlfreec.seg.gz
wget https://s3.amazonaws.com/kf-openaccess-us-east-1-prd-pbta/data/release-v2-20190809/pbta-gene-expression-kallisto.rds
wget https://s3.amazonaws.com/kf-openaccess-us-east-1-prd-pbta/data/release-v2-20190809/pbta-histologies.tsv
wget https://s3.amazonaws.com/kf-openaccess-us-east-1-prd-pbta/data/release-v2-20190809/pbta-snv-mutect2.vep.maf.gz
wget https://s3.amazonaws.com/kf-openaccess-us-east-1-prd-pbta/data/release-v2-20190809/pbta-snv-strelka2.vep.maf.gz
wget https://s3.amazonaws.com/kf-openaccess-us-east-1-prd-pbta/data/release-v2-20190809/pbta-sv-manta.tsv.gz
```

Example of download one single file
```
## aws-cli
aws s3 cp s3://kf-openaccess-us-east-1-prd-pbta/data/release-v2-20190809/pbta-histologies.tsv /path/to/local/folder/

## curl
curl -O https://s3.amazonaws.com/kf-openaccess-us-east-1-prd-pbta/data/release-v2-20190809/pbta-histologies.tsv
```
