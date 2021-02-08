# SNV hotspots

In this analysis we will evaluate snv and indels from strelka2,mutect2,vardict and lancet in [MAF format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) for hotspots and recurrence.  

- To capture recurrently mutated Oncogene and Tumor suppressor genes we will be filtering the mafs using a reference gene list created from various public sources and internal review, you can find them described [here](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/fusion_filtering#inputs-used-as-reference). 


### Steps of analysis

`00-setup_db.py` sets up a database with 4 callers are tables , this script was modified from [snv-callers](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers) module to include IMPACT, Protein_position columns in all callers. 

`01-reccurence-hotspot-overlap.Rmd` filters and combines calls from all callers and summarises the data to save recurrent mutated hotspots and additionally annotates the site as seen in MSKCC cancer hotspot database (file in input folder)

   
### Run

```
bash run_recurrence_hotspot.sh 

```

