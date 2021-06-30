## Fusion prioritization

**Module authors :** Krutika Gaonkar ([@kgaonkar6](https://github.com/kgaonkar6)), Jaclyn Taroni ([@jaclyn-taroni](https://github.com/jaclyn-taroni)), Jo Lynne Rokita ([@jharenza](https://github.com/jharenza)), and Komal Rathi ([@komalsrathi](https://github.com/komalsrathi))


This analysis is designed to filter artifacts and annotate fusion calls from STARfusion and Arriba fusion callers with the goal of prioritizing oncogenic fusions. 
We considered all inframe and frameshift fusion calls with a minimum of 1 junction reads and at least one gene partner expressed (FPKM > 1) to be potential true calls. 
We then removed fusion calls that had many spanning fragment reads compared to junction reads (`spanning fragment read count` minus `junction read count` greater than ten) as potential false positives. 
We retained fusion calls if the fused genes were detected by both callers, the same fusion was recurrent within a cancer_group (>2 samples), the fusion was specific to the cancer_group. 
We removed calls for which one gene was 5' or 3' fused to more than five different other genes within a sample as potential false positives. 
We annotated putative driver fusions and prioritized fusions when at least one fused gene was a known kinase, oncogene, tumor suppressor, curated transcription factor, on the COSMIC Cancer Gene Census list.
We also annotated fusions between pairs of genes where the fusion was observed in TCGA.
If a fusion was annotated as a putative oncogenic fusion, it was retained if it was detected by either caller.
We also gather counts for recurrent fusions and fused genes found in more than 3 participants per cancer group and represent them as binary matrices per sample.

#### Inputs from data download
* fusion-starfusion.tsv.gz : aggregated starfusion calls; a column tumor_id with the samples BS ID is added to each sample files
* fusion-arriba.tsv.gz : aggregated arriba calls; a column tumor_id with the samples BS ID is added to each sample files ; a column annots is added from running FusionAnnotator
* gene-expression-rsem-tpm-collapsed.rds

#### Inputs used as reference
* genelistreference.txt and fusionreference.txt formatted in code [here](https://gist.github.com/kgaonkar6/02b3fbcfeeddfa282a1cdf4803704794): 

Annotation | File | Source  
------ | ---------- | --------- 
| pfamID                        | http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/pfamDesc.txt.gz     | UCSC pfamID Description database |
| Domain Location               | http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ucscGenePfam.txt.gz | UCSC pfamID Description database |
| TCGA fusions                  | https://tumorfusions.org/PanCanFusV2/downloads/pancanfus.txt.gz             | TumorFusions: an integrative   resource for cancer-associated transcript fusions PMID: 29099951  |
| Transcription Factors | Table S1 https://ars.els-cdn.com/content/image/1-s2.0-S0092867418301065-mmc2.xlsx | @doi:10.1016/j.cell.2018.01.029
| Oncogenes                     | http://www.bushmanlab.org/assets/doc/allOnco_Feb2017.tsv                    | www.bushmanlab.org |
| Tumor suppressor genes (TSGs) | https://bioinfo.uth.edu/TSGene/Human_TSGs.txt?csrt=5027697123997809089      | Tumor Suppressor Gene Database   2.0 PMIDs: 23066107, 26590405 |
| Kinases                       | http://kinase.com/human/kinome/tables/Kincat_Hsap.08.02.xls |      The protein kinase complement of the human genome PMID: 12471243 |
| COSMIC genes                  | https://cancer.sanger.ac.uk/census | Catalogue of Somatic Mutations   in Cancer |
| Pediatric-specific oncogenes  | _MYBL1, SNCAIP, FOXR2, TTYH1, TERT_ | doi:10.1073/pnas.1300252110,   doi:10.1038/nature11327, doi:10.1016/j.cell.2016.01.015, doi:10.1038/ng.2849,   doi:10.1038/ng.3438, doi:10.1002/gcc.22110, doi:10.1016/j.canlet.2014.11.057,   doi:10.1007/s11910-017-0722-5 |
| Pediatric-specific TSGs | _BCOR_, _QKI_  | doi:10.1016/j.cell.2016.01.015, doi:10.1038/ng.3500 |

IGH-@,IGH@ , IGL-@ and IGL@ were also added to reference list as oncogenic genes because StarFusion output contains these gene symbols instead of IGL/IGH as per public databases.


* Brain_FPKM_hg38_matrix.txt.zip : GTex brain samples FPKM data
The code to generate genelistreference.txt and fusionreference.txt is available here: https://gist.github.com/kgaonkar6/02b3fbcfeeddfa282a1cdf4803704794#file-format_reference_gene_list-r


#### Outputs saved to data download
* fusion-putative-oncogenic.tsv
* fusion-recurrent-fusion-bycancergroup.tsv
* fusion-recurrent-fusion-bysample.tsv
* fusion-recurrently-fused-genes-bycancergroup.tsv
* fusion-recurrently-fused-genes-bysample.tsv

### Run script
use OPENPBTA_BASE_SUBTYPING=1 to run this module using the histologies-base.tsv from data folder while running molecular-subtyping modules for release.
```sh
OPENPBTA_BASE_SUBTYPING=1 run_fusion_merged.sh 
```

OR by default uses histologies.tsv from data folder
```sh
bash run_fusion_merged.sh
```


#### Order of scripts in analysis
`01-fusion-standardization.R` : Standardizes fusion calls from STARFusion and Arriba

`02-fusion-filtering.R` : Filters artifacts by removing readthroughs and NEIGHBORS from annots column; annots column also contains red flag databases that are used to further filter common/normal occuring fusions; filters all fusions where both genes are expressed FPKM < 1 are removed as non-expressed; requires fusions to have at least 1 JunctionSpanningRead and SpanningFragCount-JunctionReadCount to be <100

`03-Calc-zscore-annotate.R` : Calculates z-score for gene fused gene's expression compared to GTeX brain samples and annotates differential expression status

`04-project-specific-filtering.Rmd` : Performs project specific filtering. We prioritize the fusions as putative-oncogenic fusions if any fused gene in the fusion is annotated as kinases, oncogenes, tumor suppressors, curated transcription factors or present in COSMIC Cancer Gene Census list. We also annotated fusions if they are present in TCGA fusions list.
All fusion calls are additionally have columns `reciprocal_exists` to specify if within the Sample a fusion GeneX--GeneY has a reciprocal GeneY--GeneX . `DomainRetainedGene1A` and `DomainRetainedGene1B` are added to identify kinase domain retention for Gene1A (5` Gene) and Gene1B (3` Gene).
Oncogene annotated fusions do not need to be in both callers to be retained however if these fusions are found in more than 4 cancer groups we treat them as false calls and remove them.
To scavenge back non-oncogenic fusions that are recurrently found uniquely in a cancer group we kept fusions that were called by both callers and if >2 samples per hcancer group called the fusion.
We removed the non-oncogenic fusions with genes fused more than 5 times in a samples or found in more than 1 cancer group as potential artifact. 

`05-QC_putative_onco_fusion_dustribution.Rmd` : Plots fusions found in multiple (more than 4) cancer groups in scratch/fusion-putative-oncogenic-preQC.tsv from 04-project-specific-filtering.Rmd and removes fusion calls found in more than 4 cancer groups as QC filtering.

`06-recurrent-fusions-per-cancer-group.R` : Identifies recurrent fusions and genes that are recurrently observed in fusions. We identified RNA-seq samples that can be used independently for each patient. After the selection of samples we identify which fusions and genes are recurrent (found in >3 participants per cancer group) in our `fusion-putative-oncogenic.tsv` dataset.
