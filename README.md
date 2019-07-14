# OpenPBTA-analysis

The Open Pediatric Brain Tumor Atlas (OpenPBTA) Project is an effort to describe the landscape of tumors in the [Children's Brain tumor Tissue Consortium](https://cbttc.org/) and the PNOC003 DIPG clinical trial from the [Pediatric Pacific Neuro-oncology Consortium](http://www.pnoc.us/).
This is an open analysis effort that is organized on GitHub.
There is a [companion OpenPBTA-manuscript repository](https://github.com/AlexsLemonade/OpenPBTA-manuscript/) that is being used to author a collaborative manuscript describing the effort.
The project maintainers include scientists from [Alex's Lemonade Stand Foundation's Childhood Cancer Data Lab](https://www.ccdatalab.org/) and the [Center for Data-Driven Discovery in Biomedicine at the Children's Hospital of Philadelphia](https://d3b.center/).

## How to Participate

### Existing Planned Analyses

### Proposing a New Analysis

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

### Germ Line Variant Data (Germline SNVs)

### Somatic Single Nucleotide Variant Data (Somatic SNVs)

### Somatic Copy Number Variant Data (Somatic CNVs)

### Gene Expression Estimates

### Clinical Data

## How to Add an Analysis

### Docker Container

### Folder Structure

### Analysis Script Numbering

### Output Expectations

### Continuous Integration (CI)

### Adding Analyses to CI
