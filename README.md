# Kufareva Lab rxncon and Boolnet scripts
Various scripts for the rxncon/boolnet pipeline of the Kufareva lab.

## Installation

### Getting the repository
`git clone https://github.com/palmerito0/kboolnet.git`

### Check out the proper branch
Branches:

- `master`: Working version of the pipeline (currently a WIP)
- ~~`devel`: Main development branch, code changes should go here so they can be reviewed before merging into `master`.~~ For now, keep everything in `master`

To check out a branch, run: `git checkout <branch>`

### Install dependencies

R dependencies: `ggplot2`, `dplyr`, `openxlsx`, `googledrive`,  `tidyr`, `numbers`, `xml2`, `BoolNet`, `egg`, `optparse`, `plyr`, and `RCy3`
To install these dependencies, enter an R environment run: `install.packages(pkgs=c("ggplot2", "dplyr", "openxlsx", "googledrive",  "tidyr", "numbers", "xml2", "BoolNet", "egg", "optparse", "plyr", "RCy3"), dependencies=TRUE)`

Python dependencies: `rxncon`, `openpyxl`, `networkx`

To install these dependencies, run: `pip install rxncon openpyxl networkx`

## Contributing
A small intro to git that should be enough to get you up and running with the repository.
[See this helpful guide for more info.](https://rogerdudler.github.io/git-guide/)
Make sure to comment well!

### Making changes
After you have checked out the `devel` branch, you can make any changes you wish.
To save these changes, you must stage them, commit them, and then push them to the remote repository.

First, type `git status` to see which changes have been made to the repo and what needs to be staged. To stage a change, run: `git add <filename(s)>`

Then, run `git commit -m "<message>"` to "commit" your changes to the repository.

Finally, run `git push origin <branch>` to upload your changes to the remote repository. Only committed changes are uploaded.

### Updating your local repo
Run `git pull` to update your local repostory with all remote changes.
Git should do the job of merging changes to the same file. If you run into conflicts, shoot me an email and I'll help you resolve it. 

## Usage

[See the wiki for documentation!](https://github.com/palmerito0/kboolnet/wiki)
