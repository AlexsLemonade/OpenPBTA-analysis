# DEPRECATED: D3B TMB analysis

This *deprecated* analysis was originally used to compare these TMB calculations to what is calculated in the [`snv-callers`](../../snv-callers/README.md) module.
Its results do not end up in the final figures and results so it has since become deprecated.
It was last up to date in the v17 data release.

###### Authors : Teja Koganti for D3B

This analysis computes tumor mutation burden for different disease types.
The first script takes single MAF file and filters variants based on `Types of mutations counted` from [Friends of Cancer research](https://jitc.bmj.com/content/8/1/e000147#DC1).
TMB is computed as `(filtered variant counts* 1000000) / target BED length`.

### Calculate TMB
  1. Reads the MAF file and filters based on `Variant_Classification` column

  2. Reads in config file and defines metadata specific fields/columns and variants list etc.

  3. Reads in a target config file and determines `experimental strategy` and `cohort` from metadata file.
  MAF variants  are filtered within this target file and the BED length is calculated

  4. Calculates TMB
      - Calculates TMB based on `((# of variants)*1000000) / size of BED)`
      - Prints out the `Samplename`, `experimental_strategy`, `cohort`, `disease`, `count`, `bedlength of target BED`, `TMB` for every sample

    `Usage`: 01_calculate_tmb_targetflexible_withbothexperstrt_and_cohort.py
       [-h] -i MAF -m METADATAFILE -o OUTFILENAME -c CONFIGFILE -w
       TARGETCONFIG

        optional arguments:
          -h, --help            show this help message and exit
          -i MAF, --maf MAF     path to the MAF file
          -m METADATAFILE, --metadatafile METADATAFILE
                        path to the metadata/histology file
          -o OUTFILENAME, --outfilename OUTFILENAME
                        Out file name
          -c CONFIGFILE, --configfile CONFIGFILE
                        calculate_tmb.cfg.txt file with columns for disease,
                        samplename, variant types etc.
          -w TARGETCONFIG, --targetconfig TARGETCONFIG
                        File with experimental strategy and path to BED file

     Example -
     python3 analyses/tmb-compare/TMB_d3b_code/code/01_calculate_tmb_targetflexible.py  
        -i data/pbta-snv-consensus-mutation.maf.tsv.gz
        -m analyses/tmb-compare/TMB_d3b_code/inputs/pnoc003_and_pnoc008_and_cbttc_v17_candidate.tsv
        -o analyses/tmb-compare/TMB_d3b_code/pbta-snv-consensus.TMB.OUT.txt
        -c analyses/tmb-compare/TMB_d3b_code/config_files/calculate_tmb.cfg.json  
        -w analyses/tmb-compare/TMB_d3b_code/config_files/target_cfg.targetcombos.txt

### Plot TMB scores

    1. Takes an input file that has `Tumor_Sample_Barcode`, `cohort` and `TMB`
    2. Using matplotlib module to implement cumulative distribution function plot for every disease type
    3. Uses minimum number of samples under each disease to filter out disease types  
    4. Calculates the median line for each disease type

      `Usage` : 02_cumulative_freq_TMBplot.py [-h] -t TMB_SCORES -o OUTFILENAME -s
                                     MINSAMPLESTOPLOT

          optional arguments:
            -h, --help            show this help message and exit
            -t TMB_SCORES, --tmb_scores TMB_SCORES
                        file with TMB scores
            -o OUTFILENAME, --outfilename OUTFILENAME
                        Name of the out plot, no extension
            -s MINSAMPLESTOPLOT, --minsamplestoplot MINSAMPLESTOPLOT
                        Minimum samples from each histology/disease to plot
