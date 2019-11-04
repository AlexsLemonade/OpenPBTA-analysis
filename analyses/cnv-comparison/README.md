## Compare CNV callers

This analysis is **DEPRECATED**. 
It was designed to compare results from CNVkit and ControlFreeC when both methods produced SEG files and when the following [was noted in the README](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/0c2d0d25c01dcbbbd63f94b064a69afc9dc44ea8#data-caveats):

> We noticed ControlFreeC does not properly handle aneuploidy well for a subset of samples in that it calls the entire genome gained. 

As of `release-v7-20191031`, the CNV files are in two different formats (see: [CNVkit format](https://cnvkit.readthedocs.io/en/stable/fileformats.html) and [ControlFreeC format](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/format/controlfreec-tsv.md)).

Consensus copy number calls are tracked here: https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/128
