#! /bin/bash

set -eu

export MAKE="make -j2"

Rscript -e "install.packages( pkgs = c( 'devtools', 'optparse' ), repos = 'https://www.stats.bris.ac.uk/R/', lib = '${R_LIBS_USER}' )"
Rscript -e "install.packages( 'BiocManager', repos = 'https://www.stats.bris.ac.uk/R/', lib = '${R_LIBS_USER}' )"
Rscript -e "BiocManager::install( c( 'DNAcopy','pROC','PRROC', 'graphics'), ask = FALSE, lib = '${R_LIBS_USER}'  )"
Rscript -e "devtools::install_github( 'francescojm/CRISPRcleanR@v2.2.1', ask = FALSE, lib = '${R_LIBS_USER}', dependencies = T )"

# Need to consider how to call scripts within the pipeline when not using the Docker image
# Submodule build not installing dependencies
Rscript -e "devtools::install_github( 'cancerit/RCRISPR@develop', ask = FALSE, lib = '${R_LIBS_USER}', dependencies = T )"
