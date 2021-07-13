# Independent Samples

## Summary

Many analyses that involve mutation frequencies or co-occurence require that all samples be independent.
However, the PBTA+GMKF data set includes many cases where multiple speciments were taken from a single individual.
This analysis creates lists of samples such that there are no cases where more than one specimen is included from each individual.

As different analyses may require different sets of data, we actually generate a few different sets, stored in the `results` subdirectory. We also run the analyses based on different 'independent_level', either 'each-cohort' or 'all-cohorts'. When running with 'each-cohort', we call independent samples (based on Kids_First_Participant_ID) for each cohort+cancer_type - and same samples (based on Kids_First_Participant_ID) in different cohorts are called "independent". When running with 'all-cohorts', we call independent samples (based on Kids_First_Participant_ID) regardless of cohort or cancer_type - and same samples (based on Kids_First_Participant_ID) in different cohorts are considered the same.

The following output are generated when we run with 'all-cohorts'
* Primary specimens only with whole genome sequence (WGS):  
`independent-specimens.wgs.primary.tsv`
* Relapse specimens with WGS:  
`independent-specimens.wgs.relapse.tsv`
* Primary and relapse specimens with WGS:  
`independent-specimens.wgs.primary-plus.tsv`
* Primary specimens only with either WGS or whole exome sequence (WXS) or Panel:  
`independent-specimens.wgswxspanel.primary.tsv`
* Relapse specimens only with either WGS or whole exome sequence (WXS) or Panel:  
`independent-specimens.wgswxspanel.relapse.tsv`
* Primary and relapse specimens with WGS or WXS or Panel:  
`independent-specimens.wgswxspanel.primary-plus.tsv`

The following output are generated when we run with 'each-cohort'
* Primary specimens only with whole genome sequence (WGS):  
`independent-specimens.wgs.primary.eachcohort.tsv`
* Relapse specimens with WGS:  
`independent-specimens.wgs.relapse.eachcohort.tsv`
* Primary and relapse specimens with WGS:  
`independent-specimens.wgs.primary-plus.eachcohort.tsv`
* Primary specimens only with either WGS or whole exome sequence (WXS) or Panel:  
`independent-specimens.wgswxspanel.primary.eachcohort.tsv`
* Relapse specimens only with either WGS or whole exome sequence (WXS) or Panel:  
`independent-specimens.wgswxspanel.relapse.eachcohort.tsv`
* Primary and relapse specimens with WGS or WXS or Panel:  
`independent-specimens.wgswxspanel.primary-plus.eachcohort.tsv`

Simiarly, for independent RNA sapmles, we also run with either 'all-cohorts' or 'each-cohort'.
When run with 'each-cohort', independent DNA samples ran with 'each-cohort' was used as starting point (see code for details) and when run with 'all-cohorts', independent DNA samples ran with 'all-cohorts' was used as starting point.

The following output are generated when we run with 'all-cohorts'
* Primary and relapse RNA-Seq specimens matching WGS/WXS/Panel independent sample_ids plus only-RNA-Seq 
`independent-specimens.rnaseq.primary-plus.tsv`
`independent-specimens.rnaseq.primary.tsv`
`independent-specimens.rnaseq.relapse.tsv`

The following output are generated when we run with 'each-cohort'
* Primary and relapse RNA-Seq specimens matching WGS/WXS/Panel independent sample_ids plus only-RNA-Seq 
`independent-specimens.rnaseq.primary-plus.eachcohort.tsv`
`independent-specimens.rnaseq.primary.eachcohort.tsv`
`independent-specimens.rnaseq.relapse.eachcohort.tsv`

## Generating sample lists

To generate the independent sample lists and associated analysis of redundancies in the overall data set, run the following script from the project root directory:

use OPENPBTA_BASE_SUBTYPING=1 to run this module using the pbta-histologies-base.tsv from data folder while running molecular-subtyping modules for release.
```sh
OPENPBTA_BASE_SUBTYPING=1 ../analyses/independent-samples/run-independent-samples.sh 
```

OR by default uses histologies.tsv from data folder
```sh
bash analyses/independent-samples/run-independent-samples.sh
```

## Methods
When presented with more than one specimen from a given individual with a specific cancer group and cohort, the script randomly selects one specimen to include, with preference for primary tumors and whole genome sequences where available.
There is also a preference for the earliest collected samples, but as this data is not currently available, that code is currently deleted.

When multiple RNA-Seq samples exist per participant, the script matches the independent whole genome or whole exome sample_ids to gather matched RNA-Seq sample. If participant has onle RNA-Seq sample then a primary (and relapse if applicable) sample is randomly selected per participant per cancer group per cohort. 

## Relevant links
The methods are described in the manuscript here:
 https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#selection-of-independent-samples

 Output data files are also described in the main README here:
 https://github.com/AlexsLemonade/OpenPBTA-analysis#data-formats
