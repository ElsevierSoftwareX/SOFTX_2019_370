#!/bin/bash

#set -x
#set -e

# Instructions for building on OSX
# The OpenMS source will be installed in $HITIME_BASE/opems
# you can put it elsewhere, adjust paths accordingly
#
#
#1) Install dependencies
#
# - Xcode: https://developer.apple.com/xcode/
#
# - other dependencies using homebrew:

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

echo "brew installing dependencies"

brew install autoconf
brew install automake
brew install libtool
brew install cmake
brew install qt4
brew tap homebrew/science
brew tap homebrew/versions
#brew install libsvm xerces-c boost coinmp eigen
brew install libsvm xerces-c boost eigen

echo "building OPENMS contrib packages"

cd $HITIME_BASE
mkdir openms

cd openms
git clone https://github.com/OpenMS/contrib.git
mkdir contrib-build
cd contrib-build
cmake -D CMAKE_CXX_COMPILER=clang++ -D CMAKE_C_COMPILER=clang -D BUILD_TYPE=GLPK $HITIME_BASE/openms/contrib
cmake -D CMAKE_CXX_COMPILER=clang++ -D CMAKE_C_COMPILER=clang -D BUILD_TYPE=SEQAN $HITIME_BASE/openms/contrib
cmake -D CMAKE_CXX_COMPILER=clang++ -D CMAKE_C_COMPILER=clang -D BUILD_TYPE=WILDMAGIC $HITIME_BASE/openms/contrib
cmake -D CMAKE_CXX_COMPILER=clang++ -D CMAKE_C_COMPILER=clang -D BUILD_TYPE=GLPK $HITIME_BASE/openms/contrib
cmake -D CMAKE_CXX_COMPILER=clang++ -D CMAKE_C_COMPILER=clang -D BUILD_TYPE=EIGEN $HITIME_BASE/openms/contrib

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
cd $HITIME_BASE
git clone https://github.com/bjpop/HiTIME-CPP
#switch to the branch that you want
cd HiTIME-CPP
git checkout -b threads origin/threads

cd $HITIME_BASE/HiTIME-CPP/score

#cmake -D OpenMS_DIR=$HITIME_BASE/openms/openms_build/ .
cmake -D OpenMS_DIR=$HITIME_BASE/openms/openms_build/ -D CMAKE_PREFIX_PATH="$HITIME_BASE/openms/contrib-build;/usr;/usr/local" .
make
