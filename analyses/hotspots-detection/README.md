# Gather SNV hotspots

In this analysis we will evaluate snv and indels from strelka2,mutect2,vardict and lancet in [MAF format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) for overlaps within hotspots in MSKCC cancer hotspots v2 file downloaded [here](https://github.com/kgaonkar6/OpenPBTA-analysis/blob/recurrence-snv/analyses/hotspots-detection/input/hotspots_v2.xls) as well as use a TERT promoter region overlap to gather known mutations.

## Hotspot database for filtering
- [hotspots_v2.xls](https://www.cancerhotspots.org/files/hotspots_v2.xls) is divided into snv and indels in 2 tabs in the excel file provided in their [database](https://www.cancerhotspots.org/#/download) we will combine the 2 sheets to create a `hotspot_database_amino_acid` dataframe 

AND

- A genomic region hotspots will be in `hotspot_database_genomic` to check of TERT promoter overlaps:
TERT promoter region were gathered from  [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4852159/).
>In 2013, two hotspot point mutations were found in the TERT promoter in 71% of melanomas (32,33). The mutations were located 124bp and 146bp upstream of the translation start site and referred to as C228T and C250T, respectively, based on their hg19 genomic coordinates.

However, since we don't have the `C228T and C250T` annotation for the upstream mutations for TERT mutations, we idenitified the genomic locations in literature for hg38 genome in [this paper](https://www.mdpi.com/1422-0067/21/17/6034/htm) :
> These mutations occur at two positions upstream of the transcription starting site, at −124 bp (nucleotide polymorphism G > A, g.1295228 (chr5, 1, 295, 228 assembly GRCh37) or g.1295113 (chr5, 1, 295, 113 assembly GRCh38)) and −146 bp (nucleotide polymorphism G > A, g.1295250 (chr5, 1, 295, 250, assembly GRCh37) or g.1295135 (chr5, 1, 295, 135 assembly GRCh38)) 

- chr5 position 1295113, annotated as existing variant `rs1242535815`,`COSM1716563`,`COSM1716558`,  is 66bp away from TSS and corresponds to C228T
- chr5 position 1295135, annotated as existing variant `COSM1716559` is 88 bp away from TSS and corresponds to C250T promoter variant.


## Filtering calls that overlap hotspots

Each caller maf is filtered by:
- `IMPACT %in% c('HIGH', 'MODERATE', 'MODIFIER')` to remove any LOW mutations in the given amino acid position in hotspot database
- `Hugo_Symbol  %in% c(hotspot_database_amino_acid$Hugo_Symbol,hotspot_database_genomic$Hugo_Symbol)`
- Amino acid position overlap ( using a dataframe of hotspots with Hugo_Symbol and Amino_Acid_Position)
- TERT promoter region overlap (using a genomic region overlap filtering strategy )

Note: If a hotspot amino acid position has a silent mutation, such that the change in nucleotide creates the same amino acid in protein as Reference it is filtered out. Since these are hotspots we expect the amino acid sites to be canonical protein sequence for the gene since the annotation picks a canonical transcript first for annotation if available. 


### Steps of analysis

`00-subset-maf.R` filters and combines calls from all callers with filters for non-silent mutations in hotspot sites.
`01-create-hotspot-maf.Rmd` all calls that are not part of the consensus maf but overlap hotspots are scavenged back to a maf file along with consensus(3/3)calls. 
   
### Run

```
bash run_overlaps_hotspot.sh 

```
