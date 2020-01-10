# OpenPBTA-analysis

Pediatric brain tumors are the most common solid tumors and the leading cause of cancer-related death in children.
Our ability to understand and successfully treat these diseases is hindered by small sample sizes due to the overall rarity of unique molecular subtypes and tainted grouped analyses resulting from misclassification. 
In September of 2018, the [Children's Brain Tumor Tissue Consortium](https://cbttc.org/) released the [Pediatric Brain Tumor Atlas (PBTA)](https://cbttc.org/pediatric-brain-tumor-atlas/), a genomic dataset (whole genome sequencing, whole exome sequencing, RNA sequencing, proteomic, and clinical data) for nearly 1,000 tumors, available from the [Gabriella Miller Kids First Portal](https://kidsfirstdrc.org/). 

The Open Pediatric Brain Tumor Atlas (OpenPBTA) Project is a global open science initiative to comprehensively define the molecular landscape of tumors of 944 patients from the CBTTC and the PNOC003 DIPG clinical trial from the [Pediatric Pacific Neuro-oncology Consortium](http://www.pnoc.us/) through real-time, collaborative analyses and [collaborative manuscript writing](https://github.com/AlexsLemonade/OpenPBTA-manuscript/) on GitHub. 

The OpenPBTA operates on a pull request model to accept contributions from community participants.
The maintainers have set up continuous integration software to confirm the reproducibility of analyses within the project’s Docker container.
The collaborative manuscript is authored using [Manubot](https://manubot.org) software to provide an up-to-date public version of the manuscript. 
The project maintainers include scientists from [Alex's Lemonade Stand Foundation's Childhood Cancer Data Lab](https://www.ccdatalab.org/) and the [Center for Data-Driven Discovery in Biomedicine at the Children's Hospital of Philadelphia](https://d3b.center/).
We invite researchers to join OpenPBTA to help rigorously characterize the genomic landscape of these diseases to enable more rapid discovery of additional mechanisms contributing to the pathogenesis of pediatric brain and spinal cord tumors and overall accelerate clinical translation on behalf of patients. 

**New to the project? Please be sure to read the following documentation before contributing:**

1. Learn about the fundamental data used for this project in [**`doc/data-formats.md`**](./doc/data-formats.md) and [**`doc/data-files-description.md`**](./doc/data-files-description.md)
	+ A history of data releases can be found in [**`doc/release-notes.md`**](./doc/release-notes.md)
2. See what analyses are being performed in [**`analyses/README.md`**](./analyses/README.md) 
3. Read the remainder of this README document in full.
4. Read our contributing guidelines in [**`CONTRIBUTING.md`**](./CONTRIBUTING.md) in full.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Data Description](#data-description)
- [How to Obtain OpenPBTA Data](#how-to-obtain-openpbta-data)
  - [Data Access via Download Script](#data-access-via-download-script)
  - [Data Access via CAVATICA](#data-access-via-cavatica)
- [How to Participate](#how-to-participate)
  - [Join the Cancer Data Science Slack](#join-the-cancer-data-science-slack)
  - [Planned Analyses](#planned-analyses)
  - [Proposing a New Analysis](#proposing-a-new-analysis)
  - [Implementing an Analysis](#implementing-an-analysis)
    - [Analytical Code and Output](#analytical-code-and-output)
    - [Software Dependencies](#software-dependencies)
    - [Pull Request Model](#pull-request-model)
- [How to Add an Analysis](#how-to-add-an-analysis)
  - [Folder Structure](#folder-structure)
  - [Documenting Your Analysis](#documenting-your-analysis)
  - [Analysis Script Numbering](#analysis-script-numbering)
  - [Output Expectations](#output-expectations)
  - [Docker Image](#docker-image)
    - [Development in the Project Docker Container](#development-in-the-project-docker-container)
      - [RStudio](#rstudio)
  - [Local Development](#local-development)
    - [RStudio](#rstudio-1)
  - [Continuous Integration (CI)](#continuous-integration-ci)
    - [Working with the subset files used in CI locally](#working-with-the-subset-files-used-in-ci-locally)
    - [Adding Analyses to CI](#adding-analyses-to-ci)
    - [Adding Analyses with Multiple Steps](#adding-analyses-with-multiple-steps)
      - [1. File and merge a pull request for adding `01-filter-samples.R` to the repository.](#1-file-and-merge-a-pull-request-for-adding-01-filter-samplesr-to-the-repository)
      - [2. File and merge a pull request for adding `02-cluster-heatmap.R` to the repository.](#2-file-and-merge-a-pull-request-for-adding-02-cluster-heatmapr-to-the-repository)
      - [3. File and merge a pull request for the shell script that runs the entirety of `gene-expression-clustering`.](#3-file-and-merge-a-pull-request-for-the-shell-script-that-runs-the-entirety-of-gene-expression-clustering)
    - [Passing variables only in CI](#passing-variables-only-in-ci)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Data Description

The OpenPBTA dataset includes gene expression, fusion, as well as somatic mutation, copy number, structural and variant results in combined tsv or matrix format.

Below is a summary of biospecimens by sequencing strategy:


| Experimental Strategy | Normal | Tumor |
|-----------------------|--------|-------|
| Targeted DNA Panel | 1 | 1 |
| RNA-Seq | 0 | 1072 |
| WGS | 801 | 940 |
| WXS | 31 | 31 |


All sequencing was performed on nucleic acids extracted from fresh-frozen tissues using paired-end strategies.
The [manuscript methods section](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#data-generation) has additional details.

Below is a detailed table of [broad histologies](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#who-classification-of-disease-types) for the 1072 RNA-Seq biospecimens:


| Broad Histology | N |
|-----------------------------------------------|-----|
| Benign tumor | 38 |
| Choroid plexus tumor | 11 |
| CNS Embryonal tumor | 5 |
| CNS neuroblastoma | 5 |
| Diffuse astrocytic and oligodendroglial tumor | 225 |
| Embryonal tumor | 161 |
| Ependymal tumor | 93 |
| Germ cell tumor | 13 |
| Histiocytic tumor | 5 |
| Low-grade astrocytic tumor | 256 |
| Lymphomas | 1 |
| Meningioma | 29 |
| Mesenchymal non-meningothelial tumor | 21 |
| Metastatic secondary tumors | 7 |
| Neuronal and mixed neuronal-glial tumor | 80 |
| NOS Embryonal tumor | 17 |
| Other tumor | 6 |
| Pre-cancerous lesion | 14 |
| Tumor of cranial and paraspinal nerves | 44 |
| Tumor of pineal region | 5 |
| Tumors of sellar region | 36 |


Below is a table of number of tumor biospecimens by phase of therapy (DNA and RNA):


| Phase of Therapy | N |
|---------------------------------|------|
| Diagnosis/Initial CNS Tumor | 1507 |
| Progressive | 291 |
| Progressive Disease Post-Mortem | 11 |
| Recurrence | 128 |
| Second Malignancy | 35 |
| Unavailable | 72 |


## How to Obtain OpenPBTA Data

We are releasing this dataset on both [CAVATICA](https://cavatica.sbgenomics.com) and AWS S3.
Users performing analyses, should always refer to the symlinks in the `data/` directory and not files within the release folder, as an updated release may be produced before a publication is prepared.

**The data formats and caveats are described in more detail in [`doc/data-formats.md`](doc/data-formats.md).
For brief descriptions of the data files, see the [`data-files-description.md`](doc/data-files-description.md) file included in the download.**

Use the [data issue template](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/new?assignees=&labels=data&template=data-question.md&title=) to file issues if you have questions about or identify issues with OpenPBTA data.

### Data Access via Download Script

We have created a shell script that will download the latest release from AWS S3.
macOS users must install `md5sum` before running the download script the first time. 
This can be installed with [homebrew](https://brew.sh/) via the command `brew install coreutils` or [conda/miniconda](https://docs.conda.io/projects/conda/en/latest/) via the command `conda install -c conda-forge coreutils`.
_Note: the `download-data.sh` script now has the ability to skip downloads of unchanged files, but if you previously installed md5sum via brew you'll need to run `brew unlink md5sha1sum && brew install coreutils` first to take advantage of this new feature._

Once this has been done, run `bash download-data.sh` to acquire the latest release.
This will create symlinks in `data/` to the latest files.
It's safe to re-run `bash download-data.sh` to check that you have the most recent release of the data.
We will update the default release number whenever we produce a new release.

### Data Access via CAVATICA

For any user registered on CAVATICA, the latest release of OpenPBTA data can be accessed from the CAVATICA public projects below:
- [Pediatric Brain Tumor Atlas Open Access Data - CBTTC](https://cavatica.sbgenomics.com/u/cavatica/pbta-cbttc/)
- [Pediatric Brain Tumor Atlas Open Access Data - PNOC003](https://cavatica.sbgenomics.com/u/cavatica/pbta-pnoc003/)

Users downloading via CAVATICA should place the data files within a `data/release` folder and then create symlinks to those files within `/data`.


## How to Participate

### Join the Cancer Data Science Slack

<p><img style = "padding: 0 15px; float: left;" src = logo/slack-cancer-data-science-logo.png width = 75></p>
<p style="margin-top: 20px;"> </p>
<p><b>Have general questions or need help getting started using GitHub?</b>
You can join the <a href = "http://ccdatalab.org/slack"> Cancer Data Science Slack</a> to connect with OpenPBTA organizers, other project participants, and the broader cancer data science community.
Sign up and join the <strong>#open-pbta</strong> channel to get started!</p>

### Planned Analyses

There are certain analyses that we have planned or that others have proposed, but which nobody is currently in charge of completing.
Check the existing [issues](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues) to identify these.
If you would like to take on a planned analysis, please comment on the issue noting your interest in tackling the issue in question.
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

This section describes the general workflow for implementing analytical code, and more details are [described below](#how-to-add-an-analysis).
The first step is to identify an existing analysis or propose a new analysis, engage with the project maintainers to clarify the goals of the analysis, and then get the go ahead to move forward with the analysis.

#### Analytical Code and Output

You can perform your analyses via a script (R or Python) or via a notebook (R Markdown or Jupyter).
Your analyses should produce one or more *artifacts*.
Artifacts include both vector or high-resolution figures sufficient for inclusion in a manuscript as well as new summarizations of the data (tables, etc) that are intended for either use in subsequent analyses or distribution with the manuscript.

#### Software Dependencies

Analyses should be performed within the project's [Docker container](https://github.com/AlexsLemonade/OpenPBTA-analysis#docker-container).
We use a single monolithic container in these analyses for ease of use.
If you need software that is not included, please edit the Dockerfile to install the relevant software or file a [new issue on this repository](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/new) requesting assistance.

#### Pull Request Model

Analyses are added to this repository via [Pull Requests](https://github.com/AlexsLemonade/OpenPBTA-analysis/compare).
**Please read the [Pull Request section of the contribution guidelines](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/CONTRIBUTING.md#pull-requests) carefully.**
We are using continuous integration software applied to the supplied test datasets to confirm that the analysis can be carried out successfully within the Docker container.

## How to Add an Analysis

Users performing analyses, should **always** refer to the symlinks in the `data/` directory and not files within the release folder, as an updated release may be produced before a publication is prepared.

### Folder Structure

Our folder structure is designed to separate each analysis into its own set of notebooks that are independent of other analyses.
Within the `analyses` directory, create a folder for your analysis.
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
 We expect that the file suffix (`.csv`, `.tsv`) accurately denotes the format of the added files.
The author has also included a `README.md` ([see Documenting Your Analysis](#documenting-your-analysis)).

```
OpenPBTA-analysis
├── CONTRIBUTING.md
├── README.md
├── analyses
│   ├── existing-analysis-1
│   └── new-analysis
│       ├── 01-preprocess-data.Rmd
│       ├── 02-run-analyses.Rmd
│       ├── 03-make-figures.Rmd
│       ├── README.md
│       ├── plots
│       │   ├── figure1.pdf
│       │   └── figure2.pdf
│       ├── results
│       │   └── tabular_summary.csv
│       └── run-new-analysis.sh
├── data
└── scratch
```

### Documenting Your Analysis

A goal of the OpenPBTA project is to create a collection of workflows that are commonly used for atlas papers.
As such, documenting your analytical code via comments and including information summarizing the purpose of your analysis is important.

When you file the first pull request creating a new analysis module, add your module to the [Modules At A Glance table](analyses#modules-at-a-glance).
This table contains fields for the directory name, what input files are required, a short description, and any files that you expect other analyses will rely on.
As your analysis develops and input or output files change, please check this table remains up to date. 
This step is included in the pull request reproducibility checklist.

When an analysis module contains multiple steps or is nearing completion, add a `README.md` file that summarizes the purpose of the module, any known limitations or required updates, and includes examples for how to run the analyses to the folder.

### Analysis Script Numbering

As shown above, analysis scripts within a folder should be numbered from `01` and are intended be run in order.
If the script produces any intermediate files, these files should be placed in `../../scratch`, which is used as described above.
A shell script that runs all analytical code in the intended order should be added to the analysis directory (e.g. `run-new-analysis.sh` above).
See the [continuous integration instructions for adding analyses with multiple steps](#adding-analyses-with-multiple-steps) for more information.

### Output Expectations

The CI system that we use will generate, as artifacts, the contents of the `analyses` directory applied over a small test dataset.
Our goal is to capture all of the outputs that will be used for the [OpenPBTA-manuscript](https://github.com/AlexsLemonade/OpenPBTA-manuscript/) as artifacts.
Files that are primarily graphic should be placed in a `plots` subdirectory of the analysis's folder.
Files that are primarily tabular results files should be placed in a `results` subdirectory of the analysis's folder.
Files that are intermediate, which means that they are useful within an analysis but do not provide outputs intended for tables, figures, or supplementary tables or figures of the [OpenPBTA-manuscript](https://github.com/AlexsLemonade/OpenPBTA-manuscript/), should be placed in `../../scratch`.

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

**If you are a Mac or Windows user, the default limit for memory available to Docker is 2 GB.
You will likely need to increase this limit for local development.**
[[Mac documentation](https://docs.docker.com/docker-for-mac/#resources), [Windows documentation](https://docs.docker.com/docker-for-windows/#advanced)]

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

### Local Development

While we encourage development within the Docker container, it is also possible to conduct analysis without Docker if that is desired.
In  this case, it is important to ensure that local or personal settings such as file paths or installed packages and libraries are not assumed in the analysis.

#### RStudio

We have supplied an RStudio project (`OpenPBTA-analysis.Rproj`) file at the root of the project to aid in organization and encourage reproducible defaults for analysis.
In particular, we do not source `.Rprofile` files in new sessions or save/restore workspaces.

### Continuous Integration (CI)

We use continuous integration (CI) to ensure that the project Docker image will build if there are any changes introduced to the [`Dockerfile`](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/Dockerfile) and that all analysis code will execute.

We have put together data files specifically for the purpose of CI that contain all of the features of the full data files for only a small subset of samples.
You can see how this was done by viewing [this module](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/create-subset-files).
We use the subset files to cut down on the computational resources and time required for testing.

Provided that your analytical code points to the symlinks in the `data/` directory per [the instructions above](#how-to-add-an-analysis), adding the analysis to the CI (see below) will run your analysis on this subset of the data.
Do not hardcode sample names in your analytical code: there is no guarantee that those samples will be present in the subset files.

#### Working with the subset files used in CI locally

If you would like to work with the files used in CI locally, e.g., for debugging, you can obtain them from AWS by running the following in the root directory of the project:

```
bash scripts/download-ci-files.sh
```

Running this will change the symlinks in `data` to point to the files in `data/testing`.

#### Adding Analyses to CI

For an analysis to be run in CI, it must be added to the Circle CI configuration file, [`.circleci/config.yml`](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/.circleci/config.yml).
A new analysis should be added as the last step of the `run_analyses` section.

Here is an example analysis that simply lists the contents of the data directory that contains the files for the test:

```
      - run:
          name: List Data Directory Contents
          command: ./scripts/run_in_ci.sh ls data/testing
```

Using `./scripts/run_in_ci.sh` allows you to run your analyses in the project Docker container.

If you wanted to add running an Rscript called `cluster-samples.R` that was in an analysis folder called `gene-expression-clustering`, you would add this script to continuous integration with:

```
      - run:
          name: Cluster Samples
          command: ./scripts/run_in_ci.sh Rscript analyses/gene-expression-clustering/cluster-samples.R
```

This would run the `cluster-samples.R` on the subset files that are specifically designed to be used for CI.

#### Adding Analyses with Multiple Steps

There is a different procedure for adding an analysis comprised of multiple scripts or notebooks to CI.
Per [the contribution guidelines](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/CONTRIBUTING.md#size-and-composition-of-pull-requests), each script or notebook should be added via a separate pull request.
For each of these pull requests, the individual script or notebook should be added as its own run in the `.circleci/config.yml` file.
This validates that the code being added can be executed at the time of review.

Once all code for an analysis has been reviewed and merged, a final pull request for the analysis that is comprised of the following changes should be filed:

* A shell script that will run all script and/or notebooks in the analysis module.
* The multiple runs from the module that are in the `config.yml` file are replaced with a single run that runs the shell script.

If the `gene-expression-clustering` analysis above instead required two scripts run sequentially (`01-filter-samples.R` and `02-cluster-heatmap.R`), we would follow the procedure below.

##### 1. File and merge a pull request for adding `01-filter-samples.R` to the repository.

In this pull request, we would add the following change to `.circleci/config.yml`.

```
      - run:
          name: Filter Samples
          command: ./scripts/run_in_ci.sh Rscript analyses/gene-expression-clustering/01-filter-samples.R
```

##### 2. File and merge a pull request for adding `02-cluster-heatmap.R` to the repository.

In this pull request, we would add the following change to `.circleci/config.yml`.
This would be added _below_ the `Filter Samples` run.

```
      - run:
          name: Cluster Samples and Plot Heatmap
          command: ./scripts/run_in_ci.sh Rscript analyses/gene-expression-clustering/02-cluster-heatmap.R
```

##### 3. File and merge a pull request for the shell script that runs the entirety of `gene-expression-clustering`.

In this pull request, we would add a shell script that runs `01-filter-samples.R` and `02-cluster-heatmap.R`.
Let's call this shell script `run-gene-expression-clustering.sh` and place it in the analysis directory `analyses/gene-expression-clustering`.

The contents of `analyses/gene-expression-clustering/run-gene-expression-clustering.sh` may look like:

```
#!/bin/bash
# This script runs the gene-expression-clustering analysis
# Author's Name 2019

set -e
set -o pipefail

Rscript --vanilla analyses/gene-expression-clustering/01-filter-samples.R
Rscript --vanilla analyses/gene-expression-clustering/02-cluster-heatmap.R
```

We would remove the runs `Filter Samples` and `Cluster Samples and Plot Heatmap` from `.circleci/config.yml` and instead replace them with a single run:

```
      - run:
          name: Cluster Samples and Plot Heatmap
          command: ./scripts/run_in_ci.sh bash analyses/gene-expression-clustering/run-gene-expression-clustering.sh
```

#### Passing variables only in CI

The analyses run in CI use only a small portion of the data so that tests can be run quickly.
For some analyses, there will not be enough samples to fully test code without altering certain parameters passed to methods.
The preferred way to handle these is to run these analyses through a shell script that specifies default parameters using environment variables.
The default parameters should be the ones that are most appropriate for the full set of data.
In CI, these will be replaced.

We might decide that it makes the most sense to run an analysis using a more permissive statistical significance threshold in CI so that some "significant" pathways still appear and subsequent code that examines them can be tested.
We'd first write code capable of taking command line parameters.
In R, we could use `optparse` to specify these in a script - imagine it's called `pathway_sig.R` and it contains an option list:
```
option_list <- list(
  optparse::make_option(
    c("-a", "--alpha"),
    type = "double",
    help = "pathway significance threshold",
  )
)
```

Then we would create a shell script (perhaps `run_pathway_sig.sh`) that uses a default environment variable. If `OPENPBTA_PATHSIG` is defined, it will be used. Otherwise, a value of 0.05 is used.
Note: the `-` before the `0.05` below is necessary notation for a default parameter and *not* designating a negative 0.05.
```
PATHSIG=${OPENPBTA_PATHSIG:-0.05}

Rscript analyses/my-path/pathway_sig.R --alpha $PATHSIG
```

We can override this by passing environment variables in `.circleci/config.yml`.
For testing, we might want to use an alpha level of 0.75 so that at least some "significant" pathways appear, which allows testing subsequent code that depends on them.
The run command in the `.circleci/config.yml` is used to specify these parameters.
```
- run:
    name: run pathway significance tests
    command: OPENPBTA_PATHSIG=0.75 ./scripts/run_in_ci.sh bash analyses/my-path/run_pathway_sig.sh
```

In this example `OPENPBTA_PATHSIG=0.75` species an environment variable `OPENPBTA_PATHSIG` that is set to 0.75.
Any environment variables prefixed with `OPENPBTA_` are passed to the specified shell script.
Environment variables without this prefix are not passed.

<!--TODO: Add instructions for running scripts from anywhere in the project?-->
