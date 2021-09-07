# c-sar

[![cancerit](https://circleci.com/gh/cancerit/C-SAR.svg?style=svg)](https://circleci.com/gh/cancerit/C-SAR)

[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)

## Notes

### CRISPRcleanR v1.2.1

The version of CRISPRcleanR in the GitHub [DESCRIPTION](https://github.com/francescojm/CRISPRcleanR/blob/master/DESCRIPTION) does not match the latest release version [1.2.1](https://github.com/francescojm/CRISPRcleanR/releases/tag/v2.2.1) which is defined in [build/install_R_packages.sh](https://gitlab.internal.sanger.ac.uk/casm/team113/nextflow_pipeines/c-sar/-/blob/develop/build/install_R_packages.sh).

## Developer notes

### Version numbers

The version is in 2 files, please ensure it is updated on each release:

- `Dockerfile`
- `nextflow.config`

### Submodules

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

#### Editing the submodule

You can edit a submodule directly in this checkout and commit back to the original repository, however it is generally
less confusing to consider the submodule as a read-only component.  If you want to push changes made in this way please
see the official git documentation linked above.

#### Updating the submodule

Generally submodules are only updated to specific commits and tags.  To change which commit your branch refers to:

```bash
cd submodules/rcrispr
git pull
git checkout ${TAG_OR_COMMITREF}
cd ../../
git status
# you will see the submodule is listed as dirty
git commit -m "Updated RCRISPR submodule to ${TAG_OR_COMMITREF}" submodules/rcrispr
```

### pre-commit

Please install a pre-commit to your user environment (not a venv):

```bash
curl https://pre-commit.com/install-local.py | python3 -
# may need to add to PATH

# in repo folder, enable pre-commit hooks
pre-commit install
```

### Updating licence headers

Please use [skywalking-eyes](https://github.com/apache/skywalking-eyes).

Expected workflow:

```bash
# recent build, change to apache/skywalking-eyes:0.2.0 once released
export DOCKER_IMG=ghcr.io/apache/skywalking-eyes/license-eye
```

1. Check state before modifying `.licenserc.yaml`:
   - `docker run -it --rm -v $(pwd):/github/workspace $DOCKER_IMG header check`
   - You should get some 'valid' here, those without a header as 'invalid'
1. Modify `.licenserc.yaml`
1. Apply the changes:
   - `docker run -it --rm -v $(pwd):/github/workspace $DOCKER_IMG header fix`
1. Add/commit changes

The check is executed in the CI pipeline which will fail if expected files are missing the license.

*DO NOT* edit the header in the files, please modify the date component of `content` in `.licenserc.yaml`.  The only files needing manual update being:

- `README.Rmd`

If you need to make more extensive changes to the license carefully test the pattern is functional.

## LICENSE

```
Copyright (c) 2021 Genome Research Ltd

Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>

This file is part of C-SAR.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

1. The usage of a range of years within a copyright statement contained within
this distribution should be interpreted as being equivalent to a list of years
including the first and last year specified and all consecutive years between
them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
2009, 2011-2012’ should be interpreted as being identical to a statement that
reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
2009, 2010, 2011, 2012’.
```
