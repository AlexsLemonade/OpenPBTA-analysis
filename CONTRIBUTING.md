# Contribution guidelines for the OpenPedCan-analysis

## Issues

The best place to initially engage with the project is likely to be via creation of new issues or participation in existing issues at <https://github.com/PediatricOpenTargets/ticket-tracker>.
The OpenPedCan-analysis currently does not have a GitHub Issues tab, because it is a fork repository from the [`AlexsLemonade/OpenPBTA-analysis` repository](https://github.com/AlexsLemonade/OpenPBTA-analysis).

All developers and participants can submit issues for the following purposes:

- Propose a new analysis module.
- Update an existing analysis module.
- Address limitations of the current data processing approaches.
- Address difficulties interacting with the data.
- Address other roadblocks.
- Discuss analyses.

Participants who wish to perform an analysis of the data as part of this effort should either identify an existing, planned analysis that they wish to tackle or propose a new one.

## Pull requests

Contributions to the analysis repository operate on a pull request model.
We expect participants to actively review pull requests, with a particular focus on pull requests that include analyses within their areas of expertise.

_If you are new to git and GitHub, check out [this Hello World guide](https://guides.github.com/activities/hello-world/), [the GitHub documentation for cloning a repository](https://help.github.com/en/articles/cloning-a-repository), [the GitHub documentation for changes from a remote repository](https://help.github.com/en/articles/getting-changes-from-a-remote-repository) and [the GitHub documentation for reviewing changes in pull requests](https://help.github.com/en/articles/reviewing-changes-in-pull-requests)._

### Filing a pull request from your own branch

Here, we include a _typical_ workflow for filing a pull request to contribute to this project.

1. Fork this repository to your own account: https://help.github.com/en/articles/fork-a-repo
2. Treat `PediatricOpenTargets/OpenPedCan-analysis` as an upstream repository.
You should follow these instructions for configuring a remote that points to an upstream repository here: https://help.github.com/en/articles/configuring-a-remote-for-a-fork
3. When you want to make changes and ultimately file a pull request against the `dev` branch of this repository, you can start by making a branch from the `dev` branch of your fork named for the changes you intend to make. 
For example, you could name your branch `gene-expression-clustering` 
Then, [check out the new branch](https://gist.github.com/markSci5/5916003) (e.g., `gene-expression-clustering`).
Commit and push your changes to this branch. 
You will then be able to file a pull request from this branch: https://help.github.com/en/articles/creating-a-pull-request-from-a-fork
4. Once your pull request is approved and merged by the project maintainers, it's time to keep the `dev` branch of your fork up-to-date with the `dev` branch of this repository by following this process: https://help.github.com/en/articles/syncing-a-fork.
5. Once your pull request is merged by the project maintainer, please go to the repo, find the branch that the PR is based on, and delete the branch from the repo. 
Here's a typical series of steps for performing that:

```
# fetch the changes from the PediatricOpenTargets/OpenPedCan-analysis repository
git fetch upstream

# checkout your own dev branch
git checkout dev

# merge the upstream/dev branch into your own dev branch
git merge upstream/dev

# push the changes to your dev branch to the remote repository
git push
```

You do not want to commit changes to your `dev` branch, any changes to your `dev` branch should come from the upstream `dev` branch.

### Size and composition of pull requests

An ideal pull request is small enough for reviewers to review the code in detail and focused on a single area or, if adding a new file entirely, a single file.
Implementing an analysis will often require more than one notebook or script to be added to the repository.
It is best to submit _multiple pull requests_ for these analyses rather than a single, large pull request when an analysis is completed.
This facilitates scientific discussion and reduces the burden on reviewers (see [Peer review](#peer-review)).

The recommended size of a pull request is less than 400 lnes of code.
If your pull request contains more than 400 lines of code, it is recommended to break the pull request into multiple smaller stacked pull requests, which could be implemented with the instructions in the [Creating stacked pull requests](creating-stacked-pull-requests) section.

As the author of a pull request, consider what reviewers who have not been working on the analysis need to know to perform an effective review. 
We've put together a [pull request template](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/blob/dev/.github/PULL_REQUEST_TEMPLATE.md) to help.

#### Creating stacked pull requests

If you find that the changes on your branch include multiple files and a large number of lines of code have been changed, you may want to file a [**draft pull request**](https://help.github.com/en/articles/about-pull-requests#draft-pull-requests) to get feedback on splitting up the changes into multiple pull requests.
You can read more about how to split up pull requests [here](https://graysonkoonce.com/stacked-pull-requests-keeping-github-diffs-small/) and in the ["4 Git strategies for Pull Requests splitting" section](https://www.thedroidsonroids.com/blog/splitting-pull-request#4-git-strategies-for-pull-requests-splitting) of this article.
We include a simple example, adapted from the first link, below.

I have a new analysis on a branch called `new-analysis` that includes three scripts: `01-first-script.R`, `02-second-script.R`, and `03-third-script.R` and a bash script `run-analysis.sh` that provides the commands to run each script.

All three scripts have been committed to the `new-analysis` branch.
To file three pull requests, one for each script, I could take the following approach:

**Add the first script to its own branch `new-analysis-first`.**

We're creating the new branch with all the changes from `new-analysis` unstaged.

```sh
git checkout -b new-analysis-first new-analysis
git reset dev
```

Now we're ready to add, commit, and push `01-first-script.R` and the bash script `run-analysis.sh` with the information to run `01-first-script.R`.

```sh
# add the first script
git add analyses/new-analyses/01-first-script.R
git add analyses/new-analyses/run-analysis.sh

# stash all other changes (e.g., the second and third script)
git stash --include-untracked --keep-index

# commit + push the changes to the first script
git commit -m "Add first script for new analysis"
git push origin new-analysis-first
```

**We're then ready to file a pull request from the `new-analysis-first` branch.**

To get a pull request ready for `02-second-script.R` and an updated `run-analysis.sh` with information to run `02-second-script.R`, we'd do the following:

```sh
# create a new branch for this purpose
git checkout -b new-analysis-second

# pop the stash that contains the second and third scripts
git stash pop

# add, stash, commit, push
git add analyses/new-analyses/02-second-script.R
git add analyses/new-analyses/run-analysis.sh
git stash --include-untracked --keep-index
git commit -m "Add second script for new analysis"
git push origin new-analysis-second
```

These steps can be repeated for `03-third-script.R` and `run-analysis.sh` can be updated at each step.

### Updating your out-of-date branch

The dev branch for `PediatricOpenTargets/OpenPedCan-analysis` is protected such that pull request branches must be up-to-date before they are merged into dev.
In addition, it may be helpful to update your branch before it is approved and ready to be merged. 
For example, it is possible that a continuous integration failure is related to your branch being out-of-date (e.g., the data download script has been updated) and not due to the changes you are adding in your pull request.

You can update your branch via the GitHub interface using the `Update branch` button that appears on your pull request when this occurs. 

![](https://user-images.githubusercontent.com/6311784/38159810-cd0f73b8-34b0-11e8-928c-0ffad49b10f6.png)

**If there are merge conflicts that need to be resolved, please allow a project organizer to resolve them for you.**

If for some reason the `Update branch` button does not work, you can do the follow from the command line (provided you have [set PediatricOpenTargets/OpenPedCan-analysis as an upstream remote](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/configuring-a-remote-for-a-fork)):

```
# checkout the out-of-date branch
git checkout <my-branch-for-pull-request>
# fetch the upstream changes in the PediatricOpenTargets/OpenPedCan-analysis dev branch
git fetch upstream
# merge the PediatricOpenTargets/OpenPedCan-analysis dev branch into the out of date branch
git merge upstream/dev
```

replacing `<my-branch-for-pull-request>` with your branch name without angle brackets.

Note: You will be prompted to enter a commit message via the [text editor configured for your system](https://help.github.com/en/github/using-git/associating-text-editors-with-git).

To update your branch and pull request on GitHub, push your changes with:

```
git push
```

## Authorship

The ultimate goal of this is effort is to describe the results of our analyses in a manuscript.
We will use the [ICMJE Guidelines](http://www.icmje.org/recommendations/browse/roles-and-responsibilities/defining-the-role-of-authors-and-contributors.html).
We expect authors to contribute to the overall design of the project by participating in issues, to contribute analyses via filing pull requests that integrate analyses into this repository, to contribute to the text by contributing sections and/or revisions to sections through pull requests on the [OpenPBTA-manuscript repository](https://github.com/AlexsLemonade/OpenPBTA-manuscript/).
It is important to note that, for authorship, these should be substantial intellectual contributions.

## Peer review

All pull requests will undergo peer review. The review intends to ensure that a pull request conforms to the following basic requirements.

- The code is correct.
- The results can be reproduced in a reasonable amount of time.
- The results are expected.
- The documentation describes the purpose, methods, results, input, output, and how to run the code.
- The code follows basic style guidelines.
- The code contains comments that explain non-obvious procedures.

Participants in this project should review pull requests, which can be done using [GitHub's review interface](https://help.github.com/articles/about-pull-request-reviews/ "GitHub: about pull request reviews"). Following are the basic guidelines for reviewing pull requests:

- Note the type of review you performed: did you look over the source code, did you look over the documentation, did you run the source code, did you look at and interpret the results or a combination of these?
- Suggest modifications or, potentially, directly suggest, _but do not commit changes to_, the pull request, so that the pull request conforms to the aforementioned basic requirements.
- Explain in **detail** why the suggested modifications are **necessary**, and how the modifications should be implemented specifically, so that the developers are able to follow the suggestions without misunderstanding.
- For any suggestions/comments that are related to the lines that were changed in the particular PR, make the comment underneath the code by clicking the "+" next to the code. Click "Start a review" for the first comment, and afterwards, click "Add single comment" to keep all comments within the same review. Click "Finish your review" once completed. 
- If there is any comments/suggestions that are **outside** of the code changes this particular PR, for example, any questions about the results, or any lines of code (although not changed) that seemed puzzling, mention them as comments after you click "Finish your review". The permalinks of the referenced results or code can be obtained by the following procedure: click the "Commits" tab of the pull request, browse the repository at a specific commit that contains the referenced results or code by clicking the `<>` button of the commit row, go to the referenced file in the repository to get file or line permalinks.
- If the suggested modifications are extensive refactoring or re-designing without any change on the code behaviors, results or run-time, consider submitting a new issue for the refactoring or re-designing, and discussing the priority and specific implementations in the new issue, after the pull request under review is merged, so that the pull request under review will not be kept open for lengthy discussions and commits on extensive refactoring or re-designing that are out of the scope of the original pull request addressed issue.

Before a repository maintainer merges a pull request, there must be at least one affirmative review.
If there is any unaddressed criticism or disapproval, a repository maintainer will determine how to proceed and may wait for additional feedback.

## GitHub labels for issues and pull requests

[GitHub labels](https://docs.github.com/en/issues/using-labels-and-milestones-to-track-work/managing-labels) are used to classify issues and pull requests, so all developers can communicate commonly used messages through labels without leaving or checking comments.

Commonly used **issue** labels:

- `blocked`: the labeled issue is blocked by other issues, pull requests, etc.
- `enhancement`: the labeled issue is an enhancement of an existing analysis.
- `low priority`: The labeled issue has lower priority than issues without the label.
- `question`: the labeled issue is a question.
- `pedcbioportal`: the labeled issue is related to PedcBioPortal.
- `harmonization`: the labeled issue is related to data harmonization.

Commonly used **pull request** labels:

- `merge next`: the labeled pull request should be merged next. Pull requests without the `merge next` label should **not** be merged.
- `ready for review`: the labeled pull request is ready for review.
- `work in progress`: the labeled pull request is a work in progress, so the code, documentation, and results are subject to changes.
- `low priority`: the labeled pull request has lower priority than pull requests without the label.
- `don't merge`: the labeled pull request should not be merged.

## Merging approved pull requests

If a pull request is approved, do **not** merge it immediately. Instead, wait for other developers to coordinate with you on when to merge the pull request and whether the pull request should [be updated with the latest dev branch](#updating-your-out-of-date-branch). The coordination can reduce the workload of developers for checking whether the dev branch updates by merging other pull requests will affect their own pull requests.
