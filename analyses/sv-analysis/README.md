## Structural Variation Analysis

**Module authors :** Yang Yang ([@yangyangclover](https://github.com/yangyangclover))


This analysis is designed to 
1. structural variations analysis
2. chromothripsis analysis
3. decipher signatures of structural variations
4. to understand the potential mechanism behind tumors through structural variations analysis

#### Inputs from data download
* independent-specimens.wgs.primary-plus.tsv
* pbta-sv-manta.tsv
* pbta-cnv-consensus.seg

#### Order of scripts in analysis
`01-process-sv-file.R` : This script is for generating sv file in shatterseek-read/signature-read input format

`02-shatterseek.R` : This script is for analyzing chromothripsis, using modified ShatterSeek code