#! /bin/bash
# Copyright (c) 2021 Genome Research Ltd
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
# This file is part of C-SAR.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.
#

set -xe

if [[ -z "${TMPDIR}" ]]; then
  TMPDIR=/tmp
fi

set -u

if [ "$#" -lt "1" ] ; then
  echo "Please provide an installation path such as /opt/ICGC"
  exit 1
fi

# get path to this script
SCRIPT_PATH=`dirname $0`;
SCRIPT_PATH=`(cd $SCRIPT_PATH && pwd)`

# get the location to install to
INST_PATH=$1
mkdir -p $1
INST_PATH=`(cd $1 && pwd)`
echo $INST_PATH

# get current directory
INIT_DIR=`pwd`

CPU=`grep -c ^processor /proc/cpuinfo`
if [ $? -eq 0 ]; then
  if [ "$CPU" -gt "6" ]; then
    CPU=6
  fi
else
  CPU=1
fi
echo "Max compilation CPUs set to $CPU"

SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR/distro # don't delete the actual distro directory until the very end
mkdir -p $INST_PATH/bin
cd $SETUP_DIR

# make sure tools installed can see the install loc of libraries
set +u
export LD_LIBRARY_PATH=`echo $INST_PATH/lib:$LD_LIBRARY_PATH | perl -pe 's/:\$//;'`
export LIBRARY_PATH=`echo $INST_PATH/lib:$LIBRARY_PATH | perl -pe 's/:\$//;'`
export C_INCLUDE_PATH=`echo $INST_PATH/include:$C_INCLUDE_PATH | perl -pe 's/:\$//;'`
export PATH=`echo $INST_PATH/bin:$INST_PATH/python3/bin:$PATH | perl -pe 's/:\$//;'`
export MANPATH=`echo $INST_PATH/man:$INST_PATH/share/man:$MANPATH | perl -pe 's/:\$//;'`
export PYTHONPATH=`echo $PYTHONPATH | perl -pe 's/:\$//;'`
set -u

# Python requirements
if [ ! -e $SETUP_DIR/python_requirements.success ]; then
  pip3 install --no-cache-dir -r $SCRIPT_PATH/requirements.txt
  touch $SETUP_DIR/python_requirements.success
fi

# MAGeCK
if [ ! -e $SETUP_DIR/mageck.success ]; then
  curl -sSL --retry 10 -o mageck.tar.gz https://downloads.sourceforge.net/project/mageck/${VER_MAGECK_MAJOR}/mageck-${VER_MAGECK_MAJOR}.${VER_MAGECK_MINOR}.tar.gz
  mkdir mageck
  tar --strip-components 1 -C mageck -xzf mageck.tar.gz
  cd mageck
  python3 setup.py install --prefix=$INST_PATH/python3
  cd $SETUP_DIR
  rm -rf mageck.* mageck/*
  touch $SETUP_DIR/mageck.success
fi

# BAGEL
if [ ! -e $SETUP_DIR/bagel.success ]; then
  BAGEL_INST_DIR=$INST_PATH/bagel
  mkdir $BAGEL_INST_DIR
  git clone https://github.com/hart-lab/bagel.git ${BAGEL_INST_DIR}
  cd ${BAGEL_INST_DIR}
  git reset --hard ${VER_BAGEL}
  ln -s $BAGEL_INST_DIR/*.py $INST_PATH/bin
  cd $SETUP_DIR
  touch $SETUP_DIR/bagel.success
fi

# Nextflow
if [ ! -e $SETUP_DIR/nextflow.success ]; then
  curl -sSL https://github.com/nextflow-io/nextflow/releases/download/v${VER_NEXTFLOW}/nextflow-${VER_NEXTFLOW}-all > $INST_PATH/bin/nextflow
  chmod +x $INST_PATH/bin/nextflow
  touch $SETUP_DIR/nextflow.success
fi

echo "#!/usr/bin/env bash
nextflow run $INST_PATH/c-sar/main.nf \"\$@\"
" > $INST_PATH/bin/c-sar

chmod +x $INST_PATH/bin/c-sar
