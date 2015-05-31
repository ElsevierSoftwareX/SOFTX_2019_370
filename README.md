# HiTIME-CPP

HiTIME-CPP is a program for identifying double peaks in mass spectrometry
data implemented in C++. HiTIME-CPP is designed to be integrated with 
existing [ProteoWizard](http://proteowizard.sourceforge.net) libraries. This
implementation is based on the 
[Python implementation](https://github.com/bjpop/HiTIME). 

Slides describing the HiTIME algorithm and the development of HiTIME-CPP are
available on [Slideshare](http://goo.gl/106Yvr).

## License

HiTIME is released as open source software under the terms of the 3-Clause BSD License.
See the contents of the file called LICENSE in the top level of the source
code repository for a copy of the terms.

## Build Instructions

### Step 1: Clone this repository

```
    git clone https://github.com/bjpop/HiTIME-CPP
```

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

```
    cd src
    make
```

Running `make` in the `src` directory should be sufficent to build the
`hitime.out` binary. It may be necessary to modify the `INC` and `LIBDIR`
variables in `Makefile` if your library locations are different.

### Step 4: Test

Test data is included in this repository. Running the following command
should produce meaningful output saved in out.txt:

```
    ./hitime.out ../data/testing.mzML out.txt
```

Or launch a job on the cluster using the test slurm script:

```
    cd ../scripts
    sbatch test.slurm
```

## Usage

```
    Usage:     ./hitime [-options] [arguments]
    
    options:   -h  show this help information
               -i  ratio of doublet intensities (isotope 
                   / parent)
               -r  full width at half maximum for 
                   retention time in number of scans
               -R  retention time width boundary in 
                   standard deviations
               -p  m/z tolerance in parts per million
               -m  m/z full width at half maximum in 
                   parts per million
               -M  m/z window boundary in standard 
                   deviations
               -D  m/z difference for doublets
               -s  minimum number of data points 
                   required in each sample region
               -o  turn on full output, including zero 
                   score points
    
    arguments: mzML_file     path to mzML file
               out_file      path to output file
    
    example:   ./hitime example.mzML output.txt
```

## Documentation

The HiTIME-CPP source has been documented using 
[Doxygen](http://www.stack.nl/~dimitri/doxygen/index.html) type comments. These 
can be used to produce documentation in HTML and PDF (via Latex) formats. 

[Installation instructions](http://www.stack.nl/~dimitri/doxygen/manual/install.html). 
are available on the Doxygen website. With Doxygen installed the following 
command can be run to produce the documentation.

```
    doxygen Doxyfile
```

The result should be the creation on a `docs/` directory with `html` and `latex`
subdirectories. The HTML documentation should be viewable using any browser. To
produce PDF documentation change to the `latex` directory and run:

```
    make pdf
```

A `refman.pdf` file should be created with the HiTIME-CPP documentation. Please
note PDF documentation requires PDFLatex to be available.

## Please Note:

HiTIME-CPP is currently under active development, as such functionality can
be expected to change frequently and without notice. These instructions have
been designed for the development setup and may require significant
modification on your system.
