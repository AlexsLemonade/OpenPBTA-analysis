# OpenPBTA-analysis

The Open Pediatric Brain Tumor Atlas (OpenPBTA) Project is an effort to describe the landscape of tumors in the [Children's Brain tumor Tissue Consortium](https://cbttc.org/) and the PNOC003 DIPG clinical trial from the [Pediatric Pacific Neuro-oncology Consortium](http://www.pnoc.us/).
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
