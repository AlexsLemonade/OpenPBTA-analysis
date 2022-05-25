## Running GISTIC 2.0

**DISCLAIMER:** The `run-gistic` module **requires** the OpenPBTA project Docker container ([`ccdlopenpbta/open-pbta`](https://hub.docker.com/r/ccdlopenpbta/open-pbta)). 
The project Docker container requirement makes it distinct from other modules that are able to be run locally if you have the dependencies installed.
The code in this module uses absolute paths tied to the specific installation location on the project Docker image.

**Module author:** Jaclyn Taroni ([@jaclyn-taroni](https://github.com/jaclyn-taroni))

### Running GISTIC

Running this module via the bash script will first run `scripts/prepare_seg_for_gistic.R` and generate `pbta-cnv-consensus-gistic-only.seg.gz` file that is needed as the input for GISTIC module. 
This additional step is needed since `pbta-cnv-consensus.seg.gz` contains `copy.num=NA` for any neutral call. 
Here, `scripts/prepare_seg_for_gistic.R` takes `pbta-cnv-consensus.seg.gz` as input and uses `tumor_ploidy` from the histologies file as `copy.num`, to generate `pbta-cnv-consensus-gistic-only.seg.gz`.
`pbta-cnv-consensus-gistic-only.seg.gz` is used as input to GISTIC.

```
bash run-gistic-module.sh
```

This will also generate **incomplete** GISTIC output for short histologies with greater than 100 samples.

`run-gistic-module.sh` includes steps for generating array list files.
An array list file is (from the standalone GISTIC documentation included in the installation on the Docker container):

> The array list file is an optional file identifying the subset of samples to be used in the analysis.  It is a one column file with an optional header (array).  The sample identifiers listed in the array list file must match the sample names given in the segmentation file.  

The array list files are generated using the `scripts/generate-array-file.R` R script, which takes as arguments the histologies file (`pbta-histologies.tsv`) and a SEG file and will create files that contain biospecimen IDs that meet the following criteria: the biospecimen has the user-specified `filter_value` argument in the user-specified `filter_column` argument 2) the biospecimen ID is in the SEG file that will be used for GISTIC.

Currently, we generate array list files (see the `array_list_files` directory) filtering on the `short_histology` column of `pbta-histologies.tsv`.

Compressed GISTIC output folders are in `results`.

```
results
├── pbta-cnv-consensus-gistic.zip
├── pbta-cnv-consensus-hgat-gistic.zip
├── pbta-cnv-consensus-lgat-gistic.zip
└── pbta-cnv-consensus-medulloblastoma-gistic.zip
```

##### Why is the GISTIC output incomplete?

For runs using array list files, we see the following error at the broad analysis step:

```
Focal GISTIC completed without error
Running broad analysis...
Reconstructing genome: amp
Reconstructing genome: aod
Reconstructing genome: del
Reconstructing genome: doa
Calculating median of arm values...
arm 1: 1p 11910 markers
arm 2: 1q 10078 markers
arm 3: 2p 9007 markers
arm 4: 2q 14373 markers
arm 5: 3p 8744 markers
arm 6: 3q 9943 markers
arm 7: 4p 4730 markers
arm 8: 4q 13741 markers
arm 9: 5p 4508 markers
arm 10: 5q 12758 markers
arm 11: 6p 5339 markers
arm 12: 6q 10788 markers
arm 13: 7p 5753 markers
arm 14: 7q 9655 markers
arm 15: 8p 4113 markers
arm 16: 8q 9747 markers
arm 17: 9p 3856 markers
arm 18: 9q 7202 markers
arm 19: 10p 3752 markers
arm 20: 10q 9050 markers
arm 21: 11p 5014 markers
arm 22: 11q 7816 markers
arm 23: 12p 3255 markers
arm 24: 12q 9501 markers
arm 25: 13p 551 markers
arm 26: 13q 9501 markers
arm 27: 14p 441 markers
arm 28: 14q 8545 markers
arm 29: 15p 551 markers
arm 30: 15q 7739 markers
arm 31: 16p 3138 markers
arm 32: 16q 4252 markers
arm 33: 17p 2183 markers
arm 34: 17q 5463 markers
arm 35: 18p 1414 markers
arm 36: 18q 5839 markers
arm 37: 19p 1943 markers
arm 38: 19q 2632 markers
arm 39: 20p 2526 markers
arm 40: 20q 3372 markers
arm 41: 21p 390 markers
arm 42: 21q 3241 markers
arm 43: 22p 511 markers
arm 44: 22q 3182 markers
     1

     2

     3

     4

     5

     6

     7

     8

     9

    10

    11

    12

    13

    14

    15

    16

    17

    18

    19

    20

    21

    22

    23

    24

    25

    26

    27

    28

    29

    30

    31

    32

    33

    34

    35

    36

    37

    38

    39

    40

    41

    42

    43

    44

Error using line
Vectors must be the same lengths.

Error in gistic_broad_analysis (line 209)



Error in run_gistic20 (line 130)



Error in run_gistic2_from_seg (line 249)



Error in gp_gistic2_from_seg (line 97)



MATLAB:samelen
Warning: Objects of specgraph.scattergroup class exist - not clearing this class or any of its superclasses
Warning: Objects of scribe.legendinfo class exist - not clearing this class or any of its superclasses
Warning: Objects of scribe.legendinfochild class exist - not clearing this class or any of its superclasses
Warning: Objects of scribe.legend class exist - not clearing this class or any of its superclasses
Warning: Objects of graphics.panbehavior class exist - not clearing this class or any of its superclasses
Warning: Objects of graphics.zoombehavior class exist - not clearing this class or any of its superclasses
Warning: Objects of graphics.rotate3dbehavior class exist - not clearing this class or any of its superclasses
Warning: Objects of graphics.datacursorbehavior class exist - not clearing this class or any of its superclasses
Warning: Objects of graphics.ploteditbehavior class exist - not clearing this class or any of its superclasses
```

We have not gotten to the bottom of this as of yet.