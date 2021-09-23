## Running GISTIC 2.0

**DISCLAIMER:** The `run-gistic` module **requires** the OpenPBTA project Docker container ([`ccdlopenpbta/open-pbta`](https://hub.docker.com/r/ccdlopenpbta/open-pbta)). 
The project Docker container requirement makes it distinct from other modules that are able to be run locally if you have the dependencies installed.
The code in this module uses absolute paths tied to the specific installation location on the project Docker image.

**Module author:** Jaclyn Taroni ([@jaclyn-taroni](https://github.com/jaclyn-taroni))

### Running GISTIC on the copy number consensus SEG file

The module can be run with the following script:
```
bash run-gistic-module.sh
```

Running this module will first run `scripts/prepare_seg_for_gistic.R` and generate `cnv-consensus-gistic-only.seg.gz` file that is needed as the input for GISTIC module. 
This additional step is needed since now in `cnv-consensus.seg.gz`, any neutral call will have `copy.num=NA` (instead of `copy.num=2` before). 
This change was made since calling `copy.num=NA` can give less noisy results in the focal-cn module. 
Here, the `scripts/prepare_seg_for_gistic.R` will take that `cnv-consensus.seg.gz`, use `tumor_ploidy` as `copy.num`, generate `cnv-consensus-gistic-only.seg.gz`, and run GISTIC on that file. 

By default, running the following will run GISTIC 2.0 on `cnv-consensus-gistic-only.seg.gz`.

Compressed GISTIC output folders are in `results`.

```
results
├── cnv-consensus-gistic.zip
└── cnv-consensus-gistic-only.seg.gz
```

