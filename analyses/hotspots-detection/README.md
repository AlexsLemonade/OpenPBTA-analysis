# Gather SNV hotspots

In this analysis we will evaluate snv and indels from strelka2,mutect2,vardict and lancet in [MAF format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) for overlaps within hotspots in MSKCC cancer hotspots v2 file downloaded [here](https://github.com/kgaonkar6/OpenPBTA-analysis/blob/recurrence-snv/analyses/hotspots-detection/input/hotspots_v2.xls) as well as use a TERT promoter region overlap to gather known mutations.

## Hotspot database for filtering
- [hotspots_v2.xls](https://www.cancerhotspots.org/files/hotspots_v2.xls) is divided into snv and indels in 2 tabs in the excel file provided in their [database](https://www.cancerhotspots.org/#/download) we have subsetted the files to contain only Hugo_Symbol + Amino_Acid_Position in `hotspot_database_2017_snv.tsv` and for indels we divide into Amino_Acid_Start and Amino_Acid_End encompassing the total region where indels were identified as hotspots in `hotspot_database_2017_indel.tsv` 

AND

- A genomic region hotspots will be in `hotspot_database_genomic` to check of TERT promoter overlaps:
TERT promoter region were gathered from  [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4852159/).
>In 2013, two hotspot point mutations were found in the TERT promoter in 71% of melanomas (32,33). The mutations were located 124bp and 146bp upstream of the translation start site and referred to as C228T and C250T, respectively, based on their hg19 genomic coordinates.

However, since we don't have the `C228T and C250T` annotation for the upstream mutations for TERT mutations, we idenitified the genomic locations in literature for hg38 genome in [this paper](https://www.mdpi.com/1422-0067/21/17/6034/htm) :
> These mutations occur at two positions upstream of the transcription starting site, at −124 bp (nucleotide polymorphism G > A, g.1295228 (chr5, 1, 295, 228 assembly GRCh37) or g.1295113 (chr5, 1, 295, 113 assembly GRCh38)) and −146 bp (nucleotide polymorphism G > A, g.1295250 (chr5, 1, 295, 250, assembly GRCh37) or g.1295135 (chr5, 1, 295, 135 assembly GRCh38)) 

- chr5 position 1295113, annotated as existing variant `rs1242535815`,`COSM1716563`,`COSM1716558`,  is 66bp away from TSS and corresponds to C228T
- chr5 position 1295135, annotated as existing variant `COSM1716559` is 88 bp away from TSS and corresponds to C250T promoter variant.


## Filtering calls that overlap hotspots

 1) Each maf file is filtered using combined gene list within `hotspot_database_2017_snv.tsv` and `hotspot_database_2017_indel.tsv`

 2) Each caller maf is filtered by:
    - CANONICAL=='YES' to make sure that the annotation corresponds to the canonical transcript 
    - BIOTYPE=='protein_coding to only consider protein coding transcripts
    - IMPACT %in% c('HIGH', 'MODERATE', 'MODIFIER') to remove any LOW mutations in the given amino acid position in hotspot database
    - Amino acid position (formatted from Protein_position column in maf file) matches for SNVs in `hotspot_database_2017_snv.tsv`
    - Splice sites are matched to HGVSp_Short values in `hotspot_database_2017_snv.tsv`
    - Amino acid position (formatted from Protein_position column in maf file) within the range of indel hotspots values in `hotspot_database_2017_indel.tsv`
    - TERT promoter region overlap (using a genomic region overlap filtering strategy )


### Steps of analysis

`00-subset-maf.R` filters and combines calls from all callers with filters for non-silent mutations in hotspot sites.
`01-create-hotspot-maf.Rmd` all calls that overlap hotspots are scavenged back to a maf file excluding Vardict-only calls because Vardict uniquely calls a large number (~39 million) of very low VAF mutations as discussed [here](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/snv-callers/README.md) suggesting these could be false calls. 

   
### Run

```
bash run_overlaps_hotspot.sh 

```
