# c-sar

## Version numbers

The version is in 2 files, please ensure it is updated on each release:

* `Dockerfile`
* `nextflow.config`

## Submodules

This project has a [submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) to make it easier to build the docker
image.  The initial checkout needs to be done via recursive strategy:

```bash
# new checkout
git clone --recurse-submodules git@....git

# existing checkout
cd repo
# update local metadata (new branches etc)
git fetch
# check branch, switch if necessary
git status
# pull any changes
git pull
# populate the submodule
git submodule update --init --recursive
```

### Editing the submodule

You can edit a submodule directly in this checkout and commit back to the original repository, however it is generally
less confusing to consider the submodule as a read-only component.  If you want to push changes made in this way please
see the official git documentation linked above.

### Updating the submodule

Generally submodules are only updated to specific commits and tags.  To change which commit your branch refers to:

```bash
cd submodules/RCRISPR
git checkout ${TAG_OR_COMMITREF}
cd ../../
git status
# you will see the submodule is listed as dirty
git commit -m "Updated RCRISPR submodule to ${TAG_OR_COMMITREF}" submodules/rcrispr
```

## Notes

### CRISPRcleanR v1.2.1

The version of CRISPRcleanR in the GitHub [DESCRIPTION](https://github.com/francescojm/CRISPRcleanR/blob/master/DESCRIPTION) does not match the latest release version [1.2.1](https://github.com/francescojm/CRISPRcleanR/releases/tag/v2.2.1) which is defined in [build/install_R_packages.sh](https://gitlab.internal.sanger.ac.uk/casm/team113/nextflow_pipeines/c-sar/-/blob/develop/build/install_R_packages.sh).
