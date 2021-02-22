# SNV hotspots

In this analysis we will evaluate snv and indels from strelka2,mutect2,vardict and lancet in [MAF format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) for hotspots and recurrently mutated sites to determine whether we are missing any pediatric brain-specific hotspots not reported in MSKCC cancer hotspots v2 file downloaded [here](https://github.com/kgaonkar6/OpenPBTA-analysis/blob/recurrence-snv/analyses/hotspots-detection/input/hotspots_v2.xls)  

- To narrow the search in cancer specific genes (Oncogene and Tumor suppressor genes) we will be filtering the mafs using a reference gene list created from various public sources and internal review, you can find them described [here](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/fusion_filtering#inputs-used-as-reference). 


### Steps of analysis

`../snv-callers/scripts/01-setup_db.py` sets up a database with 4 caller maf files as tables 

`01-combine-snv.Rmd` filters and combines calls from all callers anf filters for deleterious mutations in genes of interest for cancer and brain tumor specifically (curated by @jharenza)

   
### Run

```
bash run_recurrence_hotspot.sh 

```

