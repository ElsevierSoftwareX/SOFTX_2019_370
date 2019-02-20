[![travis](https://travis-ci.org/bjpop/HiTIME-CPP.svg?branch=master)](https://travis-ci.org/bjpop/HiTIME-CPP)

# HiTIME

HiTIME is a program for identifying twin-ion signals in
[Liquid chromatography-mass spectrometry (LCMS)](https://en.wikipedia.org/wiki/Liquid_chromatography%E2%80%93mass_spectrometry) data. 
HiTIME is designed to be integrated with existing analysis pipelines, such as those 
created with [OpenMS](https://www.openms.de/). 

HiTIME filters twin-ion signals in LCMS data. This process re-weights each data point to a Z-score of how well the point matches an idealised twin-ion signal versus alternative ion signatures.

  * The `hitime` program takes an [mzML](https://en.wikipedia.org/wiki/Mass_spectrometry_data_format#mzML) file as input and produces an mzML file as output. Intentsity values in the input file which correspond to the lower mass in a twin-ion peak are retained and scored highly in the output, all other intensities are downweighted towards zero.

## License

HiTIME is released as open source software under the terms of the [3-Clause BSD License](https://opensource.org/licenses/BSD-3-Clause).
See the contents of the file called LICENSE in the top level of the source
code repository for a copy of the terms.

## Installation

HiTIME is provided as a Docker container, which can be installed like so (assuming you have docker installed on your computer):

```
docker pull bjpop/hitime
```

Docker can be installed on all modern operating systems. Please review the Docker [installation instructions](https://docs.docker.com/engine/installation/) for more information.

## Building HiTIME from source 

HiTIME can be build from source by following the instructions in the `notes` folder:

 * Mac OS X: `notes/build.osx.sh`
 * Linux: `notes/build.linux.sh`

## Builiding HiTIME using Docker

Run this command in the top-level directory of the source tree:

```
docker build -t bjpop/hitime .
```

## Example docker usage 

We provide a convenient wrapper script to run the HiTIME docker container.

Test data is included in the `data` folder within the repository. Running the following command
should produce meaningful output saved in `results.mzML`. 

```
./hitime-docker.sh -i data/testing.mzML -o results.mzML -- -d 6.0201 -r 10 -m 230
```

Note that the wrapper script accepts arguments in a slightly different format than the regular executable program: the `-i` and `-o` arguments first, then all other arguments are specified after a `--` flag.

You might see some warnings in the output which complain about the format the of input `testing.mzML` file. You can
safely ignore them. It is just OpenMS being strict about the format of the file. 

## Usage
### Command line options

```
  -h, --help            Show this help information.
  -l, --mzlower arg     Lower M/Z offset for local max window, e.g. 0.3.
                        Must be used with 'mzupper', can not be used with
                        'mzwidth'
  -u, --mzupper arg     Upper M/Z offset for local max window, e.g. 0.7. Must
                        be used with 'mzlower', can not be used with
                        'mzwidth'
  -a, --iratio arg      Ratio of doublet intensities (isotope / parent).
                        Defaults to 1.000000
  -r, --rtwidth arg     Full width at half maximum for retention time in
                        number of scans. Eg: 10
  -m, --mzwidth arg     M/Z full width at half maximum in parts per million.
                        Eg: 230. Must be used with 'mzdelta', can not be used
                        with 'mzlower' and 'mzupper'.
  -d, --mzdelta arg     M/Z delta for doublets. Eg: 6.0201. Must be used with
                        'rtwidth' and 'mzwidth', can not be used with
                        'mzlower' and 'mzupper'.
  -z, --confidence arg  Lower confidence interval to apply during scoring (In
                        standard deviations, e.g. 1.96 for a 95% CI).
                        Default: ignore confidence intervals
      --debug           Generate debugging output
      --version         Print version number and exit
  -j, --threads arg     Number of threads to use. Defaults to 1
  -c, --cache arg       Number of input spectra to retain in cache. Defaults
                        to 50
  -i, --infile arg      Input mzML file
  -o, --outfile arg     Output mzML file
```

### Twin-ion scoring example

Note that the examples below illustrate how to call the executable program directly. If you wish to use the Docker wrapper script`hitime-docker.sh` then you will need to adjust the syntax accordingly; see the notes above about using the wrapper script. 

```
hitime -j 4 -i data/testing.mzML -o results.mzML -d 6.0201 -r 10 -m 230
```

for a computation using 4 threads, where `data/testing.mzML` contains the input mass spectrometry data in mzML format, and the output file is called `results.mzML`.

The parameters defining the taget twin-ion signal are, `-d 6.0201` the M/Z diference between the natural and heavy isotope versions of the precursor, `-r 10` the retention time (RT) full width half maximum (FWHM) size in number of RT steps (scans), `-m 230` the M/Z FWHM size in parts per million (ppm).  These values can be determined by measurement of the precursor signal in standard visulisation software.

### Local maxima example
HITIME can also be used to filter the data to only output the data point that has the largest value in a region defined by the Retention Time (RT) full width half maximum (FWHM) size, and the M/Z upper and lower bounds.  E.g.:

```
hitime -i results.mzML -o max.results.mzML -r 10 -l 0.3 -u 0.7
```
This will produce two files, `max.results.mzML` and `max.results.csv`.  The CSV file is a comma separated text file listing the local maxima.  This list can be sorted to help identify the strongest twin-ion signal matches.  The fields are RT, M/Z, score.

## Indexing your input mzML file:

HITIME assumes that the input mzML file is indexed.  To index an input file, the OpenMS `FileConverter` tool can be used, eg:

```
FileConverter -ini openms.ini -in input.mzML -out input_indexed.mzML
```

An example openms.ini can be found in `notes/openms.ini`.

## Data Sets

### Small dataset (1.9 MB)
```
curl -sL https://melbourne.figshare.com/ndownloader/files/14131079 > BSA_SUBSET.mzML
```

### Larger dataset (82 MB)

A pure protein that was reduced, reacted with heavy/light paracetamol and then digested. It contains lots of twin ions. Because these twin-ions are mostly peptide derivatives, they come in multiple charge states which effectively gives us a few different doublet spacings. There are many separated by 6.0201, a good number separated by 3.01005 and a few separated by 2.0067. We have removed any points with intensity <1000 to make the file smaller.

```
curl -sL https://melbourne.figshare.com/ndownloader/files/14131085 > BSA_FULL.mzML
```
