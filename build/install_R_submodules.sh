#! /bin/bash

set -eu

Rscript -e "install.packages('./submodules/rcrispr', repos = NULL, type = 'source', lib = '${R_LIBS_USER}', dependencies = TRUE)"
