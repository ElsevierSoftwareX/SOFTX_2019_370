# HiTIME-CPP

HiTIME-CPP is a program for identifying double peaks in mass spectrometry
data implemented in C++. HiTIME-CPP is designed to be integrated with 
existing [OpenMS](http://open-ms.sourceforge.net/) libraries. 

## License

HiTIME is released as open source software under the terms of the 3-Clause BSD License.
See the contents of the file called LICENSE in the top level of the source
code repository for a copy of the terms.

## Build Instructions

#### Step 1: Build OpenMS from source

- [Linux](http://ftp.mi.fu-berlin.de/pub/OpenMS/release-documentation/html/install_linux.html)
- [OS X](http://ftp.mi.fu-berlin.de/pub/OpenMS/release-documentation/html/install_mac.html)

See below for more notes on building OpenMS.

#### Step 2: Clone this repository

```
cd $HOME/code
git clone https://github.com/bjpop/HiTIME-CPP
```

#### Step 3: Configure the build with cmake

```
cd $HOME/code/HiTIME-CPP/src/
cmake -D OpenMS_DIR=/absolute/path/to/OpenMS-build/ .
```
**NOTE** in the command above the text `/absolute/path/to/OpenMS-build/` should be the full (absolute) path
to the place where the directory `OpenMS-build` was created. In our instructions below it would be the expansion
of `$HOME/code/openms/OpenMS-build`. For some strange reason it seems that Cmake does not like it if you have
shell variables, such as `$HOME` in the `OpenMS_DIR` variable setting.

### Step 3: Compile

```
make
```

### Step 4: Test

Test data is included in this repository. Running the following command
should produce meaningful output saved in out.txt:

```
cd data
../src/hitime -i testing.mzML -o results.mzML
```

## Usage

```
HiTIME allowed options:
  -h [ --help ]         Show this help information
  -a [ --iratio ] arg   Ratio of doublet intensities (isotope / parent)
  -r [ --rtwidth ] arg  Full width at half maximum for retention time in number
                        of scans
  -t [ --rtwindow ] arg Retention time width boundary in standard deviations
  -p [ --ppm ] arg      M/Z tolerance in parts per million
  -m [ --mzwidth ] arg  M/Z full width at half maximum in parts per million
  -z [ --mzwindow ] arg M/Z window boundary in standard deviations
  -d [ --mzdelta ] arg  M/Z delta for doublets
  -n [ --mindata ] arg  Minimum number of data points required in each sample 
                        region
  --debug               Minimum number of data points required in each sample 
                        region
  -j [ --threads ] arg  Number of threads to use
  -i [ --infile ] arg   Input mzML file
  -o [ --outfile ] arg  Output mzML file
```

For example:

```
hitime -j 4 -i example.mzML -o results.mzML 
```

for a computation using 4 threads, where `example.mzML` contains the input mass spectrometry data in mzML format, and the output file is called `results.mzML`

## Building OpenMS

These instructions assume you are building the code in `$HOME/code`. If you wish to build it somewhere else replace 
`$HOME/code` with your own location.

OpenMS requires the following build tools and libraries:

 - autoconf
 - automake
 - libtool
 - cmake
 - Qt

You will also need the gcc C compiler installed, and an internet connection (the build process downloads files).

#### Optionally specify where the gcc and g++ compilers are

On some systems `gcc` and `g++` are not in standard paths. Cmake seems to have difficulty with this.
This issue can be resolved by exporting the `CC` and `CXX` environment variables:

```
export CC=`which gcc`
export CXX=`which g++`
```

### Linux build 
 
On Linux systems the build tools and libraries can be installed via the package management system, for example on Ubuntu:

#### Install packaged dependencies

```
sudo apt-get install autoconf automake libtool cmake qt-sdk 
```

#### Create a directory to build OpenMS

```
mkdir $HOME/code/openms
cd $HOME/code/openms
```

#### Clone OpenMS contrib libraries (bundled dependencies)

```
git clone https://github.com/OpenMS/contrib.git
```

#### Make a directory to build the contrib libraries

```
mkdir $HOME/code/openms/contrib-build
cd $HOME/code/openms/contrib-build
```

#### build the ZLIB bundled contrib separately.

We have to do this on its own due to a bug in the build software. See [this ticket](https://github.com/OpenMS/contrib/issues/5).
You ought to be able to build with `-DBUILD_TYPE=ALL` directly, but it fails because of the above bug. 

```
cmake -DBUILD_TYPE=ZLIB ../contrib
```

#### build all the remaining contrib libraries

 - SEQAN, LIBSVM, XERCESC, BOOST, GSL, COINOR, BZIP2, GLPK, EIGEN, WILDMAGIC

this will take a while

```
cmake -DBUILD_TYPE=ALL ../contrib
```

#### clone the OpenMS repository

```
cd $HOME/code/openms
git clone https://github.com/OpenMS/OpenMS.git
```

#### create a build directory for OpenMS

```
mkdir $HOME/code/openms/OpenMS-build
cd $HOME/code/openms/OpenMS-build
```

#### configure the OpenMS build with cmake


```
cmake -DCMAKE_PREFIX_PATH="$HOME/code/openms/contrib-build" ../OpenMS
```

If the configuration works successfully you should see a lot of output which contains the message "You have successfully configured OpenMS." near the end.


#### build OpenMS


```
make
```

This will take a while. When it is done you should have a large number of executables in the directory `$HOME/code/openms/OpenMS-build`, such as TOPPView. 

It will also try to build the documentation. This can fail if you do not have a graphics display working (such as X11). If building the documentation fails then you can run:

```
make doc_minimal
```

### OS X build

On OS X they can be install via Homebrew:

```
brew install autoconf automake libtool cmake qt4
```

You will also need [Xcode](https://developer.apple.com/xcode/)
