## Fusion prioritization

**Module authors :** Krutika Gaonkar ([@kgaonkar6](https://github.com/kgaonkar6)), Jaclyn Taroni ([@jaclyn-taroni](https://github.com/jaclyn-taroni)), Jo Lynne Harenza ([@jharenza](https://github.com/jharenza)), and Komal Rathi ([@komalsrathi](https://github.com/komalsrathi))


This analysis is designed to filter artifacts and annotate fusion calls from STARfusion and Arriba fusion callers with the goal of prioritizing oncogenic fusions. 
We considered all inframe and frameshift fusion calls with a minimum of 1 junction reads and at least one gene partner expressed (FPKM > 1) to be potential true calls. 
We then removed fusion calls that had many spanning fragment reads compared to junction reads (`spanning fragment read count` minus `junction read count` greater than ten) as potential false positives. 
We retained fusion calls if the fused genes were detected by both callers, the same fusion was recurrent within a broad_histology (>2 samples), the fusion was specific to the broad_histology. 
We removed calls for which one gene was 5' or 3' fused to more than five different other genes within a sample as potential false positives. 
We annotated putative driver fusions and prioritized fusions when at least one fused gene was a known kinase, oncogene, tumor suppressor, curated transcription factor, on the COSMIC Cancer Gene Census list.
We also annotated fusions between pairs of genes where the fusion was observed in TCGA.
If a fusion was annotated as a putative oncogenic fusion, it was retained if it was detected by either caller.
We also gather counts for recurrent fusions and fused genes found in more than 3 participants per histology and represent them as binary matrices per sample.

#### Inputs from data download
* pbta-fusion-starfusion.tsv.gz : aggregated starfusion calls; a column tumor_id with the samples BS ID is added to each sample files
* pbta-fusion-arriba.tsv.gz : aggregated arriba calls; a column tumor_id with the samples BS ID is added to each sample files ; a column annots is added from running FusionAnnotator
* pbta-gene-expression-rsem-fpkm.polya.rds : aggregated polya samples fpkm data
* pbta-gene-expression-rsem-fpkm.stranded.rds : aggregated stranded fpm data

#### Inputs used as reference
* genelistreference.txt : known kinases, oncogenes, tumor suppressors, curated transcription factors [@doi:10.1016/j.cell.2018.01.029], COSMIC Cancer Gene Census list[https://cancer.sanger.ac.uk/census] . MYBL1 [@doi:10.1073/pnas.1300252110], SNCAIP [@doi:10.1038/nature11327], FOXR2 [@doi:10.1016/j.cell.2016.01.015], TTYH1 [@doi:10.1038/ng.2849], and TERT [@doi:10.1038/ng.3438; @doi:10.1002/gcc.22110; @doi:10.1016/j.canlet.2014.11.057; @doi:10.1007/s11910-017-0722-5] were added to the oncogene list and BCOR [@doi:10.1016/j.cell.2016.01.015] and QKI [@doi:10.1038/ng.3500] were added to the tumor suppressor gene list based on pediatric cancer literature review
* fusionreference.txt :  known TCGA fusions
* Brain_FPKM_hg38_matrix.txt.zip : GTex brain samples FPKM data

#### Outputs saved to data download
* pbta-fusion-putative-oncogenic.tsv
* pbta-fusion-recurrently-fused-genes-byhistology.tsv
* pbta-fusion-recurrently-fused-genes-bysample.tsv

### Run script
`bash run_fusion_merged.sh` 

#### Order of scripts in analysis
`01-fusion-standardization.R` : Standardizes fusion calls from STARFusion and Arriba

`02-fusion-filtering.R` : Filters artifacts by removing readthroughs and NEIGHBORS from annots column; annots column also contains red flag databases that are used to further filter common/normal occuring fusions; filters all fusions where both genes are expressed FPKM < 1 are removed as non-expressed; requires fusions to have at least 1 JunctionSpanningRead and SpanningFragCount-JunctionReadCount to be <10

`03-Calc-zscore-annotate.R` : Calculates z-score for gene fused gene's expression compared to GTeX brain samples and annotates differential expression status

`04-project-specific-filtering.Rmd` : Performs project specific filtering. We removed fusions with genes fused more than 5 times in a samples as potential artifact. We kept fusions that were called by both callers and if >2 samples per histology called the fusion. We then prioritize the fusions as putative-oncogenic fusions if any fused gene in the fusion is annotated as kinases, oncogenes, tumor suppressors, curated transcription factors or present in COSMIC Cancer Gene Census list. We also annotated fusions if they are present in TCGA fusions list. 
Annotated fusions do not need to be in both callers to be retained.

`05-recurrent-fusions-per-histology.R` : Identifies recurrent fusions and genes that are recurrently observed in fusions. We identified RNA-seq samples that can be used independently for each patient. After the selection of samples we identify which fusions and genes are recurrent (found in >3 participants per histology) in our `pbta-fusion-putative-oncogenic.tsv` dataset.
