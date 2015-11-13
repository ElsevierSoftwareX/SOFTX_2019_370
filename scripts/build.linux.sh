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
HITME_BUILD_DIR="$BASEDIR/hitme_build"
OPEN_MS_DIR="$BASEDIR/openms"
OPEN_MS_CONTRIB_BUILD_DIR="$OPEN_MS_DIR/contrib-build"
OPEN_MS_BUILD_DIR="$OPEN_MS_DIR/OpenMS-build"
HITIME_SRC_DIR="HITME_BUILD_DIR/HiTIME-CPP/src"

echo "making hitime_build directory"
mkdir $HITME_BUILD_DIR
cd $HITME_BUILD_DIR

echo "making openms directory"
mkdir $OPEN_MS_DIR
cd $OPEN_MS_DIR

echo "cloning openms contrib repository https://github.com/OpenMS/contrib.git"
git clone https://github.com/OpenMS/contrib.git > OpenMS_contrib.git_clone.out 2>&1

echo "making openms contrib build directory"
mkdir $OPEN_MS_CONTRIB_BUILD_DIR
cd $OPEN_MS_CONTRIB_BUILD_DIR

echo "building openms contrib dependencies"
cmake -DBUILD_TYPE=ZLIB ../contrib > OpenMS.contrib_build.out 2>&1
cmake -DBUILD_TYPE=ALL ../contrib >> OpenMS.contrib_build.out 2>&1

echo "cloning openms repository git clone https://github.com/OpenMS/OpenMS.git"
cd $OPEN_MS_DIR
git clone https://github.com/OpenMS/OpenMS.git > OpenMS.git_clone.out 2>&1

echo "making build directory for openms"
mkdir $OPEMS_BUILD_DIR
cd $OPEMS_BUILD_DIR 

echo "configuring openms build"
cmake -DCMAKE_PREFIX_PATH="$OPEN_MS_CONTRIB_BUILD_DIR" ../OpenMS > OpenMS_config.out 2>&1

echo "building openms"
make > openms_build.out 2>&1

cd $HITME_BUILD_DIR
echo "cloning HiTIME-CPP repository"
git clone https://github.com/bjpop/HiTIME-CPP > HiTIME.git_clone.out 2>&1

echo "configuring hitime build" 

cd $HITIME_SRC_DIR
cmake -D OpenMS_DIR=$OPEN_MS_BUILD_DIR . > HiTIME.config.out > 2>&1

echo "building hitime"
make > HiTIME.build.out > 2>&1

echo "done"
