## Counting git contributions to a branch

### Purpose

The purpose of this module is to count the git contributions to this repo in two ways:

1. Count the number of analysis modules – directories in `analyses/` that are included in the manuscript, specifially – that individual contributors have authored commits in. 
  You can find that information in `results/module_contribution_counts.tsv` and a _less summarized_ version that lists all contributors to a module in `results/module_contributors.tsv`.
2. Count the number of commits to the entire repository individual contributors have authored.

### Running the module

To count git contributions on :warning: your current branch :warning:, you can run the following:

```sh
bash run-count-contributions.sh
```

Which runs the following steps:

* `01-count-contributions.sh` which is a shell script that generates intermediate TXT files (in `scratch/`) with the output of `git shortlog` and `git log`.
* `02-format-contributions.Rmd` which turns the output of `git shortlog` and `git log` into readable files. This notebook additionally collapses people that have their contributions divided across multiple names by the logs and is responsible for removing the analysis modules that are not included in the manuscript. See the notebook for more details!

### Checking out the correct branch

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

Now you're ready to run the module!