# HiTIME-CPP

HiTIME-CPP is a program for identifying double peaks in mass spectrometry
data implemented in C++. HiTIME-CPP is designed to be integrated with 
existing [ProteoWizard](proteowizard.sourceforge.net) libraries. This
implementation is based on the 
[Python implementation](https://github.com/bjpop/HiTIME). 

# Build Instructions

### Step 1: Clone this repository

`git clone git@github.com:lazappi/HiTIME-CPP.git`

### Step 2: Load required modules

HiTIME-CPP requires both the ProteoWizard and [Boost](www.boost.org) 
libraries as well as taking advantage newer compiler features. For a 
successful build these must be available.

```
    module load gcc/4.9.1
    module load pwiz-gcc/3.0.7069
    module load boost-gcc/1.57.0
```

### Step 3: Linker path

In order to run the compiled program the ProteoWizard libraries must be able
to be located. This is done by adding them to the `LD_LIBRARY_PATH` 
environment variable.

`export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/pwiz/3.0.7069/lib/`

### Step 4: Compile

Running `make` in the HiTIME-CPP directory should be sufficent to build the
`hitime.out` binary. It may be necessary to modify the `INC` and `LIBDIR`
variables in `makefile` if your library locations are different.

### Step 5: Test

Test data is included in this repository. Running the following command
should produce meaningful output.

`./hitime.out data/testing.mzML`

## Please Note:

HiTIME-CPP is currently under active development, as such functionality can
be expected to change frequently and without notice. These instructions have
been designed for the developement setup and may require significant
modification on your system.

