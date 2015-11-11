# HiTIME-CPP

HiTIME-CPP is a program for identifying double peaks in mass spectrometry
data implemented in C++. HiTIME-CPP is designed to be integrated with 
existing [OpenMS](http://open-ms.sourceforge.net/) libraries. 

## License

HiTIME is released as open source software under the terms of the 3-Clause BSD License.
See the contents of the file called LICENSE in the top level of the source
code repository for a copy of the terms.

## Build Instructions

### Step 1: Build OpenMS from source

- [Linux](http://ftp.mi.fu-berlin.de/pub/OpenMS/release-documentation/html/install_linux.html)
- [OS X](http://ftp.mi.fu-berlin.de/pub/OpenMS/release-documentation/html/install_mac.html)

See below for more notes on building OpenMS.

### Step 2: Clone this repository

```
    git clone https://github.com/bjpop/HiTIME-CPP
```

### Step 3: Compile

Running `make` in the `src` directory should be sufficent to build the
`hitime` binary. It may be necessary to modify the `INC` and `LIBDIR`
variables in `Makefile` if your library locations are different.

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

### Linux depdencies
 
On Linux systems they can be installed via the package management system, for example on Ubuntu:

```
sudo apt-get install autoconf automake libtool cmake qt-sdk 
mkdir $HOME/code/openms
cd $HOME/code/openms

```

You will also need the gcc C compiler installed.

### OS X dependencies

On OS X they can be install via Homebrew:

```
brew install autoconf automake libtool cmake qt4
```

You will also need [Xcode](https://developer.apple.com/xcode/)
