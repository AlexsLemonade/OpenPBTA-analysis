# SNV hotspots

In this analysis we will evaluate snv and indels from strelka2,mutect2,vardict and lancet in [MAF format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) for overlaps within hotspots in MSKCC cancer hotspots v2 file downloaded [here](https://github.com/kgaonkar6/OpenPBTA-analysis/blob/recurrence-snv/analyses/hotspots-detection/input/hotspots_v2.xls) as well as use a TERT promoter region overlap to gather known mutations.


### Steps of analysis

`../snv-callers/scripts/01-setup_db.py` sets up a database with 4 caller maf files as tables 

`01-combine-snv.Rmd` filters and combines calls from all callers with filters for non-silent mutations in hotspot sites.

   
### Run

```
bash run_overlaps_hotspot.sh 

```

