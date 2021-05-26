## Chromothripsis Analysis (in progress)

**Module authors :**
Yang Yang ([@yangyangclover](https://github.com/yangyangclover))
Modifications in progress by Laura Egolf ([@LauraEgolf](https://github.com/LauraEgolf/))

This module runs ShatterSeek, classifies chromothripsis regions, and visualizes the results.

#### Inputs from data download
* independent-specimens.wgs.primary-plus.tsv
* pbta-sv-manta.tsv
* pbta-cnv-consensus.seg

#### Order of scripts in analysis
`01-process-sv-file.R` : This script reformats SV files for ShatterSeek input.

`02-shatterseek.R` : This script runs ShatterSeek and classifies chromothripsis regions as high or low confidence based on criteria recommended by the authors. [Modifications ongoing]

More scripts will be added.