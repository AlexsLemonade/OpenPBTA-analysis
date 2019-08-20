# OpenPBTA-analysis

The Open Pediatric Brain Tumor Atlas (OpenPBTA) Project is an effort to describe the landscape of tumors in the [Children's Brain Tumor Tissue Consortium](https://cbttc.org/) and the PNOC003 DIPG clinical trial from the [Pediatric Pacific Neuro-oncology Consortium](http://www.pnoc.us/).
This is an open analysis effort that is organized on GitHub.
There is a [companion OpenPBTA-manuscript repository](https://github.com/AlexsLemonade/OpenPBTA-manuscript/) that is being used to author a collaborative manuscript describing the effort.
The project maintainers include scientists from [Alex's Lemonade Stand Foundation's Childhood Cancer Data Lab](https://www.ccdatalab.org/) and the [Center for Data-Driven Discovery in Biomedicine at the Children's Hospital of Philadelphia](https://d3b.center/).

## How to Participate

### Planned Analyses

There are certain analyses that we have planned or that others have proposed, but which nobody is currently in charge of completing.
Check the existing [issues](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues) to identify these.
We have tagged a [subset of these with the label "good first issue"](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22).
These "good first issues" are ones that have a limited scope, that can be done with the software already present on the Docker container, and that we expect to have few dependencies with other issues.
If you would like to take on these or any other existing planned analysis, please comment on the issue noting your interest in tackling the issue in question.
Ask clarifying questions to understand the current scope and goals.
Then propose a potential solution.
If the solution aligns with the goals, we will ask you to go ahead and start to implement the solution.
You should provide updates to your progress in the issue.
When you file a pull request with your solution, you should note that it closes the issue in question.

### Proposing a New Analysis

In addition to the planned analyses, we welcome contributors who wish to propose their own analyses of this dataset as part of the OpenPBTA project.
Check the existing [issues](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues) before proposing an analysis to see if something similar is already planned.
If there is not a similar planned analysis, create a new issue.
The ideal issue will describe the scientific goals of the analysis, the planned methods to address the scientific goals, the input data that is required for the planned methods, and a proposed timeline for the analysis.
Project maintainers will interact on the issue to clarify any questions or raise any potential concerns.

### Implementing an Analysis

This section describes the general workflow for implementing analytical code, and more details are [described below](https://github.com/AlexsLemonade/OpenPBTA-analysis#how-to-add-an-analysis).
The first step is to identify an existing analysis or propose a new analysis, engage with the project maintainers to clarify the goals of the analysis, and then get the go ahead to move forward with the analysis.
Analyses should be performed within the project's [Docker container](https://github.com/AlexsLemonade/OpenPBTA-analysis#docker-container).
We use a single monolithic container in these analyses for ease of use.
If you need software that is not included, please edit the Dockerfile to install the relevant software or file a [new issue on this repository](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/new) requesting assistance.
You can perform your analyses via a script (R or Python) or via a notebook (R Markdown or Jupyter).
Your analyses should produce one or more *artifacts*.
Artifacts include both vector or high-resolution figures sufficient for inclusion in a manuscript as well as new summarizations of the data (tables, etc) that are intended for either use in subsequent analyses or distribution with the manuscript.
You should file a [Pull Request](https://github.com/AlexsLemonade/OpenPBTA-analysis/compare) to contribute analyses to this repository.
We are using continuous integration software applied to the supplied test datasets to confirm that the analysis can be carried out successfully within the Docker container.

## How to Obtain OpenPBTA Data

The OpenPBTA dataset includes somatic mutation and gene expression results in combined tsv or matrix format.
We are releasing this dataset on both [CAVATICA](https://cavatica.sbgenomics.com) and AWS S3.
Users performing analyses, should always refer to the symlinks in the `data/` directory and not files within the release folder, as an updated release may be produced before a publication is prepared.

### Data Access via Download Script

We have created a shell script that will download the latest release from AWS S3.
OS X users must use [homebrew](https://brew.sh/) to install `md5sha1sum` via the command `brew install md5sha1sum` before running the download script the first time.
Once this has been done, run `bash download-data.sh` to acquire the latest release.
This will create symlinks in `data/` to the latest files.
It's safe to re-run `bash download-data.sh` to check that you have the most recent release of the data.
We will update the default release number whenever we produce a new release.

### Data Access via CAVATICA

For any user registered on CAVATICA, the latest release of OpenPBTA data can be accessed from the CAVATICA public projects below:
- [Pediatric Brain Tumor Atlas Open Access Data - CBTTC](https://cavatica.sbgenomics.com/u/cavatica/pbta-cbttc/)
- [Pediatric Brain Tumor Atlas Open Access Data - PNOC003](https://cavatica.sbgenomics.com/u/cavatica/pbta-pnoc003/)

Users downloading via CAVATICA should place the data files within a `data/release` folder and then create symlinks to those files within `/data`.

## Data Formats

The release notes for each release are provided in the `release-notes.md` file that accompanies the data files.
Somatic Single Nucleotide Variant (SNV) data are provided in [Annotated MAF format](doc/format/vep-maf.md) files for each of the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#somatic-single-nucleotide-variant-calling).
Somatic Copy Number Variant (CNV) data are provided in the [SEG format](https://software.broadinstitute.org/software/igv/SEG) for each of the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#somatic-copy-number-variant-calling).
Somatic Structural Variant Data (Somatic SV) are provided in the [Annotated Manta TSV](doc/format/manta-tsv-header.md) and [LUMPY TSV](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/format/lumpy-tsv-header.md) formats produced by the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#somatic-structural-variant-calling).
Gene expression estimates from the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#gene-expression-abundance-estimation) are provided as a gene by sample matrix.
Gene Fusions produced by the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#rna-fusion-calling-and-prioritization) are provided as [Arriba TSV](doc/format/arriba-tsv-header.md) and [STARFusion TSV](doc/format/starfusion-tsv-header.md) respectively.
[Harmonized clinical data](https://alexslemonade.github.io/OpenPBTA-manuscript/#clinical-data-harmonization) are released as tab separated values.

### Data Caveats
The clinical manifest will be updated and versioned as molecular subgroups are identified based on genomic analyses.
We noticed ControlFreeC does not properly handle aneuploidy well for a subset of samples in that it calls the entire genome gained. We recommend using these results with caution while we continue benchmarking this algorithm.


## How to Add an Analysis

Users performing analyses, should **always** refer to the symlinks in the `data/` directory and not files within the release folder, as an updated release may be produced before a publication is prepared.

### Docker Image

We build our project Docker image from a versioned [`tidyverse`](https://hub.docker.com/r/rocker/tidyverse) image from the [Rocker Project](https://www.rocker-project.org/) (v3.6.0).

To add dependencies that are required for your analysis to the project Docker image, you must alter the project [`Dockerfile`](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/Dockerfile).
Note: R packages installed on this image will be installed from an [MRAN snapshot](https://mran.microsoft.com/documents/rro/reproducibility#reproducibility) corresponding to the last day that R 3.6.0 was the most recent release ([ref](https://hub.docker.com/r/rocker/tidyverse)).

If you need assistance adding a dependency to the Dockerfile, [file a new issue on this repository](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/new) to request help.

#### Development in the Project Docker Container

The most recent version of the project Docker image, which is pushed to Docker Hub after a pull request gets merged into the master branch, can be obtained via the command line with:

```
docker pull ccdlopenpbta/open-pbta:latest
```

##### RStudio

Using `rocker/tidyverse:3.6.0` as our base image allows for development via RStudio in the project Docker container. 
If you'd like to develop in this manner, you may do so by running the following and changing `<password>` to a password of you choosing at the command line:

```
docker run -e PASSWORD=<password> -p 8787:8787 ccdlopenpbta/open-pbta:latest
```

You can change the volume that the Docker container points to either via the [Kitematic GUI](https://docs.docker.com/kitematic/userguide/) or the [`--volume` flag](https://docs.docker.com/storage/volumes/) to `docker run`.

Once you've set the volume, you can navigate to `localhost:8787` in your browser if you are a Linux or Mac OS X user. 
The username will for login will be `rstudio` and the password will be whatever password you set with the `docker run` command above.

If you are a new user, you may find [these instructions](https://github.com/AlexsLemonade/RNA-Seq-Exercises/blob/master/docker-pull.md) for a setting up a different Docker container or [this guide](https://www.andrewheiss.com/blog/2017/04/27/super-basic-practical-guide-to-docker-and-rstudio/) from Andrew Heiss helpful.

### Folder Structure

Our folder structure is designed to separate each analysis into its own set of notebooks that are independent of other analyses.
Within the analyses directory, create a folder for your analysis.
Choose a name that is unique from other analyses and somewhat detailed.
For example, instead of `gene-expression`, choose `gene-expression-clustering` if you are clustering samples by their gene expression values.
You should assume that any data files are in the `../../data` directory and that their file names match what the `download-data.sh` script produces.
These files should be read in at their relative path, so that we can re-run analyses if the underlying data change.
Files that are primarily graphic should be placed in a `plots` subdirectory.
Files that are primarily tabular results files should be placed in a `results` subdirectory.
Intermediate files that are useful within the processing steps but that do not represent final results should be placed in `../../scratch/`.
It is safe to assume that files placed in `../../scratch` will be available to all analyses within the same folder.
It is not safe to assume that files placed in `../../scratch` will be available from analyses in a different folder.

An example highlighting a `new-analysis` directory is shown below.
The directory is placed alongside existing analyses within the `analyses` directory.
In this case, the author of the analysis has run their workflows in R Markdown notebooks.
This is denoted with the `.Rmd` suffix.
However, the author could have used Jupyter notebooks, R scripts, or another scriptable solution.
The author has produced their output figures as `.pdf` files.
We have a preference for vector graphics as PDF files, though other forms of vector graphics are also appropriate.
The results folder contains a tabular summary as a comma separated values file.
We expect that the file suffix (`.csv`, `.tsv`) should accurately denote the format of the files added.

```
OpenPBTA-analysis
├── CONTRIBUTING.md
├── README.md
├── analyses
│   ├── existing-analysis-1
│   └── new-analysis
│       ├── 01.preprocess-data.Rmd
│       ├── 02.run-analyses.Rmd
│       ├── 03.make-figures.Rmd
│       ├── plots
│       │   ├── figure1.pdf
│       │   └── figure2.pdf
│       └── results
│           └── tabular_summary.csv
├── data
└── scratch
```

### Analysis Script Numbering

As shown above, analysis scripts within a folder should be numbered from `01` and are intended be run in order.
If the script produces any intermediate files, these files should be placed in `../../scratch`, which is used as described above.

### Output Expectations

The CI system that we use will generate, as artifacts, the contents of the `analyses` directory applied over a small test dataset.
Our goal is to capture all of the outputs that will be used for the [OpenPBTA-manuscript](https://github.com/AlexsLemonade/OpenPBTA-manuscript/) as artifacts.
Files that are primarily graphic should be placed in a `plots` subdirectory of the analysis's folder.
Files that are primarily tabular results files should be placed in a `results` subdirectory of the analysis's folder.
Files that are intermediate, which means that they are useful within an analysis but do not provide outputs intended for tables, figures, or supplementary tables or figures of the [OpenPBTA-manuscript](https://github.com/AlexsLemonade/OpenPBTA-manuscript/), should be placed in `../../scratch`.

### Continuous Integration (CI)

We continuous integration (CI) to ensure that all analysis code will execute and that the project Docker image will build if there are any changes introduced to the [`Dockerfile`](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/Dockerfile).

We have put together files specifically for the purpose of CI that contain the full features in a data file, but only for a small subset of samples. 
(You can see how this was done by viewing [this notebook](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/create-subset-files/01-create_subset_files.nb.html).)
We use the subset files to cut down on the computational resources and time required for testing.
Provided that your analytical code points to the symlinks in the `data/` directory per [the instructions above](#how-to-add-an-analysis), adding the analysis to the CI will run on these subset files (see below).
Note that if you hardcode sample names in your analytical code and those samples are not present in the subset files, this may cause errors.

#### Adding Analyses to CI

For an analysis to be run in CI, it must be added to the Circle CI configuration file, [`.circleci/config.yml`](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/.circleci/config.yml). 
A new analysis should be added as the last step of the `run_analyses` job.

Here is an example analysis that simply lists the contents of the data directory that contains the files for the test:

```
      - run:
          name: List Data Directory Contents
          command: ./scripts/run_in_ci.sh ls data/testing
```

Using `./scripts/run_in_ci.sh` allows you to run your analyses in the project Docker container.

If you wanted to add running an Rscript called `cluster-samples.R` that was in an analysis folder called `gene-expression-structure`, you would add this script to continuous integration with:

```
      - run:
          name: Cluster Samples
          command: ./scripts/run_in_ci.sh Rscript analyses/gene-expression-structure/cluster-samples.R
```

This would run the `cluster-samples.R` on the subset files that are specifically designed to be used for CI.
