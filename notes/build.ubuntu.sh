#!/bin/bash

#set -x
#set -e

# Instructions for building on Ubuntu 
# The OpenMS source will be installed in $HITIME_BASE/openms
# you can put it elsewhere, adjust paths accordingly
#
#
#1) Install dependencies
#
# - using the Ubuntu package manager

# You will need to define the HITIME_BASE shell variable to be
# the filepath of the place where you want to make the installation.
# 
# Most likely this should be in the same (parent) directory where you
# have a copy of the HiTIME-CPP repository.
#
# For example:
#
# HITIME_BASE=$HOME/code

HITIME_BASE=$HOME/code

if [ -z ${HITIME_BASE+x} ]; then
    echo "You must set the variable HITIME_BASE in the build shell script"
    exit 1
fi

echo "installing dependencies"

sudo apt-get install cmake g++ autoconf qt4-dev-tools patch libtool make git \
     libqt4-core libqt4-dev libqt4-gui libqt4-opengl-dev automake libqtwebkit-dev
sudo apt-get install libboost-regex-dev libboost-iostreams-dev libboost-date-time-dev libboost-math-dev \
     libsvm-dev libglpk-dev libzip-dev zlib1g-dev libxerces-c-dev libbz2-dev

echo "building OPENMS contrib packages"

cd $HITIME_BASE
mkdir openms

cd openms
git clone https://github.com/OpenMS/contrib.git
mkdir contrib-build
cd contrib-build
cmake -DBUILD_TYPE=SEQAN ../contrib
cmake -DBUILD_TYPE=WILDMAGIC ../contrib
cmake -DBUILD_TYPE=EIGEN ../contrib

echo "configuring OPENMS"

cd $HITIME_BASE/openms
git clone https://github.com/OpenMS/OpenMS.git
mkdir openms_build
cd openms_build
cmake -DCMAKE_PREFIX_PATH="$HITIME_BASE/openms/contrib-build;/usr;/usr/local" ../OpenMS

echo "building OPENMS"

make
make test

# 2) Now you can build HiTIME-CPP 
# This assumes you have already checked out the repository for HiTIME-CPP like so:
# cd $HITIME_BASE
# git clone https://github.com/bjpop/HiTIME-CPP
# switch to the branch that you want
# cd HiTIME-CPP
# git checkout -b threads origin/threads

cd $HITIME_BASE/HiTIME-CPP/src

cmake -D OpenMS_DIR=$HITIME_BASE/openms/openms_build/ -D CMAKE_PREFIX_PATH="$HITIME_BASE/openms/contrib-build;/usr;/usr/local" .
make
