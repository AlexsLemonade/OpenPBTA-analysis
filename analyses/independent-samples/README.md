# Independent Samples

## Summary

Many analyses that involve mutation frequencies or co-occurence require that all samples be independent.
However, the PBTA data set includes many cases where multiple speciments were taken from a single individual.
This analysis creates lists of samples such that there are no cases where more than one specimen is included from each individual.

As different analyses may require different sets of data, we actually generate a few different sets, stored in the `results` subdirectory:
* Primary specimens only with whole genome sequence (WGS):  
`independent-specimens.wgs.primary.tsv`
* Primary and secondary specimens with WGS:  
`independent-specimens.wgs.primary-plus.tsv`
* Primary specimens only with either WGS or whole exome sequence (WXS):  
`independent-specimens.wgswxs.primary.tsv`
* Primary and secondary specimens with WGS or WXS:  
`independent-specimens.wgswxs.primary-plus.tsv`

## Generating sample lists

To generate the independent sample lists and associated analysis of redundancies in the overall data set, run the following script from the project root directory:

```sh
bash analyses/independent-samples/run-independent-samples.sh
```


## Methods

When presented with more than one specimen from a given individual, the script randomly selects one specimen to include, with preference for primary tumors and whole genome sequences where available.
There is also a preference for the earliest collected samples, but as this data is not currently available, that code is not currently relevant.

## Relevant links
The methods are described in the manuscript here:
 https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#selection-of-independent-samples

 Output data files are also described in the main README here:
 https://github.com/AlexsLemonade/OpenPBTA-analysis#data-formats
