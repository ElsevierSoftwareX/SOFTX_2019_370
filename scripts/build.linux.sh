#!/bin/bash

# Script to build HiTIME-CPP on linux.
# Automates building of OPEN-MS as well.
#
# Run this script in the directory where you want
# to build HiTIME-CPP. It will create a sub-directory
# called hitime_build.
#
#
# You need to have these things installed:
#
# autoconf
# automake
# libtool
# cmake
# Qt
#
# At VLSCI (on barcoo) you can make them available with this command
#
# module load autoconf automake libtool-gcc cmake/3.0.2 qt-gcc/4.8.5 
#
#
# On other Linux distributions you can make the available with the
# package manager, for example on Ubuntu/Debian:
#
# sudo apt-get install autoconf automake libtool cmake qt-sdk

export CC=`which gcc`
export CXX=`which g++`

BASEDIR=`pwd`
HITME_BUILD_DIR="$BASEDIR/hitime_build"

OPEN_MS_CONTRIB_BUILD_DIR="$HITME_BUILD_DIR/contrib-build"
OPEN_MS_BUILD_DIR="$HITME_BUILD_DIR/OpenMS-build"
HITIME_SRC_DIR="$HITME_BUILD_DIR/HiTIME-CPP/src"

echo "making hitime_build directory"
mkdir $HITME_BUILD_DIR
cd $HITME_BUILD_DIR

echo "cloning openms contrib repository https://github.com/OpenMS/contrib.git"
git clone https://github.com/OpenMS/contrib.git 

echo "making openms contrib build directory"
mkdir $OPEN_MS_CONTRIB_BUILD_DIR
cd $OPEN_MS_CONTRIB_BUILD_DIR

echo "building openms contrib dependencies"
cmake -DBUILD_TYPE=ZLIB ../contrib 
cmake -DBUILD_TYPE=ALL ../contrib 

echo "cloning openms repository git clone https://github.com/OpenMS/OpenMS.git"
cd $HITME_BUILD_DIR
git clone https://github.com/OpenMS/OpenMS.git 

echo "making build directory for openms"
mkdir $OPEN_MS_BUILD_DIR
cd $OPEN_MS_BUILD_DIR

echo "configuring openms build"
cmake -DCMAKE_PREFIX_PATH="$OPEN_MS_CONTRIB_BUILD_DIR" ../OpenMS 

echo "building openms"
make doc_minimal

cd $HITME_BUILD_DIR
echo "cloning HiTIME-CPP repository"
git clone https://github.com/bjpop/HiTIME-CPP 

echo "configuring hitime build" 

cd $HITIME_SRC_DIR
cmake -D OpenMS_DIR=$OPEN_MS_BUILD_DIR . 

echo "building hitime"
make 

echo "done"
