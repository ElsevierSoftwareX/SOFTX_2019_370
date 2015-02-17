# HiTIME-CPP

HiTIME-CPP is a program for identifying double peaks in mass spectrometry
data implemented in C++. HiTIME-CPP is designed to be integrated with 
existing [ProteoWizard](https://proteowizard.sourceforge.net) libraries. This
implementation is based on the 
[Python implementation](https://github.com/bjpop/HiTIME). 

Slides describing the HiTIME algorithm and the development of HiTIME-CPP are
available on [Slideshare](http://goo.gl/106Yvr).

# Build Instructions

### Step 1: Clone this repository

`git clone git@github.com:lazappi/HiTIME-CPP.git`

### Step 2: Load required modules

HiTIME-CPP requires both the ProteoWizard and [Boost](www.boost.org) 
libraries as well as taking advantage newer compiler features. For a 
successful build these must be available. **NOTE:** HiTIME-CPP may not build
correctly with other versions of gcc such as 4.8.1.

```
    module load gcc/4.9.1
    module load pwiz-gcc/3.0.7069
    module load boost-gcc/1.57.0
```

### Step 3: Compile

Running `make` in the HiTIME-CPP directory should be sufficent to build the
`hitime.out` binary. It may be necessary to modify the `INC` and `LIBDIR`
variables in `makefile` if your library locations are different.

### Step 4: Test

Test data is included in this repository. Running the following command
should produce meaningful output saved in out.txt:

`./hitime.out data/testing.mzML out.txt`

# Documentation

The HiTIME-CPP source has been documented using 
[Doxygen](http://www.stack.nl/~dimitri/doxygen/index.html) type comments. These 
can be used to produce documentation in HTML and PDF (via Latex) formats. 

[Installation instructions](http://www.stack.nl/~dimitri/doxygen/manual/install.html). 
are available on the Doxygen website. With Doxygen installed the following 
command can be run to produce the documentation.

`doxygen Doxyfile`.

The result should be the creation on a `docs/` directory with `html` and `latex`
subdirectories. The HTML documentation should be viewable using any browser. To
produce PDF documentation change to the `latex` directory and run:

`make pdf`

A `refman.pdf` file should be created with the HiTIME-CPP documentation. Please
note PDF documentation requires PDFLatex to be available.

# Please Note:

HiTIME-CPP is currently under active development, as such functionality can
be expected to change frequently and without notice. These instructions have
been designed for the development setup and may require significant
modification on your system.
