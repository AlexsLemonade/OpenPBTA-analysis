# Molecular Subtyping of Low Grade Gliomaa

**Module authors:** Krutika Gaonkar ([@kgaonkar6](https://github.com/kgaonkar6)) Jaclyn Taroni ([@jaclyn-taroni](https://github.com/jaclyn-taroni))

In this analysis we subtype LGAT samples according to the presence/absence of BRAF fusion/BRAF V600E point mutation. 

### Preprocessing
The files in the `lgat-subset` were generated using 00-subset-files-for-LGAT.R which uses files from data release folder   

### Inputs from data download
* pbta-histologies.tsv: is used to subset samples where short_histology == LGAT 
* pbta-snv-consensus-mutation.maf.tsv.gz: is used to find samples where Hugo_Symbol==BRAF & HGVSp_Short=V600E
* pbta-fusion-recurrently-fused-genes-bysample.tsv: is used to find samples with BRAF fusions
* pbta-fusion-putative-oncogenic.tsv: is used to find samples with KIAA1549--BRAF fusions


### Run script
`bash analyses/molecular-subtyping-LGAT/run_subtyping.sh`

#### Order of scripts in analysis
00-subset-files-for-LGAT.R: generates subset of wgs LGAT sample annotated with BRAF_V600E column denoting presence="Yes" and absence="No" 
01-make-lgat-final-table.Rmd: generates final table for LGAT subtyping from metadata,fusion and consenses mutation data release files 
