# Independent Samples

## Summary

Many analyses that involve mutation frequencies or co-occurence require that all samples be independent.
However, the PBTA data set includes many cases where multiple speciments were taken from a single individual.
This analysis creates lists of samples such that there are no cases where more than one specimen is included from each individual.

As different analyses may require different sets of data, we actually generate a few different sets, stored in the `results` subdirectory:
* Primary specimens only with whole genome sequence (WGS):  
`independent-specimens.wgs.primary.tsv`
* Secondary specimens with WGS:  
`independent-specimens.wgs.secondary.tsv`
* Primary and secondary specimens with WGS:  
`independent-specimens.wgs.primary-plus.tsv`
* Primary specimens only with either WGS or whole exome sequence (WXS) or Panel:  
`independent-specimens.wgswxspanel.primary.tsv`
* Secondary specimens only with either WGS or whole exome sequence (WXS) or Panel:  
`independent-specimens.wgswxspanel.secondary.tsv`
* Primary and secondary specimens with WGS or WXS or Panel:  
`independent-specimens.wgswxspanel.primary-plus.tsv`

* Primary and secondary RNA-Seq specimens matching WGS/WXS/Panel independent sample_ids plus only-RNA-Seq 
`independent-specimens.rnaseq.primary-plus.tsv`


## Generating sample lists

To generate the independent sample lists and associated analysis of redundancies in the overall data set, run the following script from the project root directory:

use OPENPBTA_BASE_SUBTYPING=1 to run this module using the pbta-histologies-base.tsv from data folder while running molecular-subtyping modules for release.
```sh
OPENPBTA_BASE_SUBTYPING=1 ../analyses/independent-samples/run-independent-samples.sh 
```

OR by default uses pbta-histologies.tsv from data folder
```sh
bash analyses/independent-samples/run-independent-samples.sh
```

## Methods

When presented with more than one specimen from a given individual, the script randomly selects one specimen to include, with preference for primary tumors and whole genome sequences where available.
There is also a preference for the earliest collected samples, but as this data is not currently available, that code is not currently relevant.

When multiple RNA-Seq samples exist per participant, the script matches the independent whole genome or whole exome sample_ids to gather matched RNA-Seq sample. If participant has onle RNA-Seq sample then a primary (and secondary if applicable) sample is randomly selected per participant  

## Relevant links
The methods are described in the manuscript here:
 https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#selection-of-independent-samples

 Output data files are also described in the main README here:
 https://github.com/AlexsLemonade/OpenPBTA-analysis#data-formats
