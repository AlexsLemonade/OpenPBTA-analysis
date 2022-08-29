## Counting git contributions to a branch

### Purpose

The purpose of this module is to count the git contributions to this repo in two ways:

1. Count the number of analysis modules – directories in `analyses/` that are included in the manuscript, specifically – that individual contributors have authored commits in.
  You can find that information in `results/module_contribution_counts.tsv` and a _less summarized_ version that lists all contributors to a module in `results/module_contributors.tsv`.
2. Count the number of commits to the entire repository individual contributors have authored.

### Steps

The shell script `run-count-contributions.sh` runs the following steps:

* `01-count-contributions.sh` which is a shell script that generates intermediate TXT files (in `scratch/count-contributions`) with the output of `git shortlog` and `git log`.
* `02-format-contributions.Rmd` which turns the output of `git shortlog` and `git log` into readable files.
This notebook additionally collapses people that have their contributions divided across multiple names by the logs and is responsible for removing the analysis modules that are not included in the manuscript.
See the notebook for more details!
* `03-set-authorship-order.Rmd` which downloads and reads in the current manuscript metadata file (i.e., from `AlexsLemonade/OpenPBTA-manuscript` `master` branch) and sets the author order based on the following criteria (quoting from the notebook itself):

> * The first, last 4 authors, and consortia authors positions are pinned to reflect scholarship and contributor roles that are not captured in the Git history of this repository (e.g., substantive analytical code review, conceptualization, supervision, or project administration).
> * Manuscript authors that have contributed to the code base for the analysis repository -- specifically the number of analysis modules that were included in the paper that a manuscript author has contributed to -- are then ordered from the second position on in decreasing order.
> In the case of ties, the order is randomly selected.
> * Manuscript authors that did not directly contribute to the code base are then randomly ordered.
> We set a seed directly before that shuffling step, using the year as the seed, to keep a consistent order in future runs if and when the Git contributions change.

The updated metadata YAML file is ignored by this repository, but can be found at `results/metadata.yaml`.


* `04-get-author-information.Rmd` which extracts relevant author information (name, email, affiliation, and ORCID) into a TSV for use during manuscript.
It creates a file `author_information.tsv`, which is ignored by this repository but will be included in the `AlexsLemonade/OpenPBTA-manuscript` repository in `submission_info/` once the GitHub Actions workflow is run.

### GitHub Actions workflow

To keep these stats reasonably up-to-date, we use a GitHub Actions (GHA) workflow ([`.github/workflows/count-git-contributions.yml`](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/.github/workflows/count-git-contributions.yml)) to rerun the module every Wednesday at 14:00 UTC.

The workflow runs `bash run-count-contributions.sh` on the `master` branch and files a new pull request with the updated tables and HTML file from the notebook.
This PR is intended to be reviewed by an organizer for accuracy before merging.

If the author order changes – i.e., the YAML output from `03-set-authorship-order.Rmd` is different from the `master` branch of `AlexsLemonade/OpenPBTA-manuscript` – a pull request will be filed in the manuscript repository.
The updated author TSV file `submission_info/author_information.tsv` will be filed as part of this pull request.

The _scheduled_ GHA workflow is the main way we expect the module to be used, but you can also make use of the [`workflow_dispatch`](https://docs.github.com/en/actions/using-workflows/events-that-trigger-workflows#workflow_dispatch)  trigger [in the GitHub browser interface](https://docs.github.com/en/actions/managing-workflow-runs/manually-running-a-workflow).

### Running the module manually

To count git contributions on :warning: your current branch :warning:, you can run the following:

```sh
bash run-count-contributions.sh
```

#### Checking out the correct branch

**You most likely want to count contributions to the `AlexsLemonade` `master` branch!**
We include some instructions for doing this locally below.

`AlexsLemonade` should be set as an `upstream` remote, which you can check with:

```sh
git remote -v
```

Which should include the following if `AlexsLemonade` is set as `upstream`:

```sh
upstream	https://github.com/AlexsLemonade/OpenPBTA-analysis.git (fetch)
upstream	https://github.com/AlexsLemonade/OpenPBTA-analysis.git (push)
```

You likely have your own `master` branch locally, so to checkout `upstream/master` with a new local name (using `upstream-master` below), you can use the following:

```sh
git checkout -b upstream-master upstream/master
```

If you already have a local `upstream-master` branch, check it out and make sure it's up to date!

Now you're ready to run the module!

#### Committing the changes from a manual run

Once you've run `bash run-count-contributions.sh` to generate stats for the `AlexsLemonade/master` branch, you'll need to _create a new branch_ and commit the modified files to that new branch.

