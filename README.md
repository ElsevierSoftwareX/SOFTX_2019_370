# HiTIME-CPP

HiTIME-CPP is a program for identifying twin-ion signals in
[Liquid chromatography-mass spectrometry (LCMS)](https://en.wikipedia.org/wiki/Liquid_chromatography%E2%80%93mass_spectrometry) data. 
HiTIME-CPP is designed to be integrated with existing analysis pipelines, such as those 
created with [OpenMS](https://www.openms.de/). 

HiTIME-CPP filters twin-ion signals in LCMS data. Each data point in the input is scored according to the goodness-of-fit of its neighbourhood to twin Gaussians of expected dimensions and spacings. This process re-weights each data point by its likelihood of being part of a true twin ion signal. 

HiTIME-CPP bundles two programs:
  * `hitime score`: takes an [mzML](https://en.wikipedia.org/wiki/Mass_spectrometry_data_format#mzML) file as input and produces an mzML file as output. Intentsity values in the input file which correspond to the lower mass in a twin-ion peak are retained and scored highly in the output, all other intensities are downweighted towards zero.
  * `hitime max`: takes the mzML output from `hitime score` and detects local maxima intensities, these should occur at the center of the lower mass in each twin ion pair. 

## License

HiTIME is released as open source software under the terms of the [3-Clause BSD License](https://opensource.org/licenses/BSD-3-Clause).
See the contents of the file called LICENSE in the top level of the source
code repository for a copy of the terms.

## Installation

HiTIME-CPP is provided as a Docker container, which can be installed like so (assuming you have docker installed on your computer):

```
docker pull bjpop/hitime
```

Docker can be installed on all modern operating systems. Please review the Docker [installation instructions](https://docs.docker.com/engine/installation/) for more information.

## Building HiTIME-CPP from source 

HiTIME-CPP can be build from source by following the instructions in the `notes` folder:

 * Mac OS X: `notes/build.osx.sh`
 * Linux: `notes/build.linux.sh`

## Builiding HiTIME-CPP using Docker

Run this command in the top-level directory of the source tree:

```
docker build -t bjpop/hitime .
```

## Example docker usage 

We provide a convenient wrapper script to run the HiTIME-CPP docker container.

Test data is included in the `data` folder within the repository. Running the following command
should produce meaningful output saved in `results.mzML`. 

```
./hitime-docker.sh score -i data/testing.mzML -o data/results.mzML
```

You might see some warnings in the output which complain about the format the of input `testing.mzML` file. You can
safely ignore them. It is just OpenMS being strict about the format of the file. 

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

