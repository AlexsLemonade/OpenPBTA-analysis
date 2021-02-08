# SNV hotspots

This analysis we evaluate snv and indels from strelka2,mutect2,vardict and lancet in [MAF format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) to look for Oncogene and TumroSuppressorGenes recurrently mutated in 1 or more callers. This reference gene list used here was created from fusion_filtering module before and described here. 

Steps of analysis

`00-setup_db.py` sets up a database with 4 callers are tables , this script was modified from snv-callers module to include IMPACT, Protein_position columns in all callers. 

`01-reccurence-hotspot-overlap.Rmd` filters and combines calls from all callers and summarises the data to save recurrent mutated hotspots and additionally annotates the site as seen in MSKCC cancer hotspot database

   
Run

```
bash run_recurrence_hotspot.sh 

```

