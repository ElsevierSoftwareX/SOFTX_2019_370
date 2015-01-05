#!/bin/env python

# Standard libraries
import argparse

# Third Party libraries
import pymzml
import numpy as np

# Parse command line arguaments

parser = argparse.ArgumentParser(description='Convert an mzML file to NumPy arrays.')

parser.add_argument('inputFile',
                    help='mzML input data file',
                    type=str)

parser.add_argument('-o',
                    '--output',
                    help='directory for output of NumPy arrays',
                    type=str,
                    default='')

opts = parser.parse_args()

# Find the number of spectra in the file and the number of points in the
# biggest spectra in order to determine the size of the NumPy matrices

# Exclude the last spectrum as it has large values inconsistent with the
# rest of the data
# This may just be a problem with the test data we are using

ms_run = pymzml.run.Reader(opts.inputFile)
mz_lengths = []
spectrum_count = 0

for spectrum in ms_run:
    spectrum_count += 1
    mz_lengths.append(len(spectrum.mz))

max_mz_length = max(mz_lengths[:-1])

# Create empty matrices using the determined dimensions. Missing values
# will remain as NaN

mzs = np.full((max_mz_length, spectrum_count-1), np.nan, dtype='float64')
intensities = np.full((max_mz_length, spectrum_count-1), np.nan,
                      dtype='float64')
times = np.full((max_mz_length, spectrum_count-1), np.nan, dtype='float64')

# Iterate over the spectra in the file and insert values into appropriate
# locations in the matrices

ms_run = pymzml.run.Reader(opts.inputFile)
spectrum_idx = 0

while spectrum_idx < spectrum_count-1:
    spectrum = ms_run.next()
    mzs[0:len(spectrum.mz), spectrum_idx] = spectrum.mz
    intensities[0:len(spectrum.i), spectrum_idx] = spectrum.i
    times[0:len(spectrum.i), spectrum_idx] = np.array(
        spectrum['MS:1000016'] * len(spectrum.i))
    spectrum_idx += 1

# Output the NumPy arrays

file_prefix = opts.inputFile.split('/')[-1].split('.')[0]

np.save(opts.output + '/' + file_prefix + '_mz.npy', mzs)
np.save(opts.output + '/' + file_prefix + '_time.npy', times)
np.save(opts.output + '/' + file_prefix + '_intensity.npy', intensities)
