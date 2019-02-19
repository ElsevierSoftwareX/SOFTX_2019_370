#ifndef HITIME_CONSTANTS_H
#define HITIME_CONSTANTS_H

#include <math.h>

//! Convert standard deviation to FWHM
const float std_dev_in_fwhm = 2.355;
//! Default MZ boundary sigma.
const float default_mz_sigma = 2;
//! Default ratio of peak intensities.
const float default_intensity_ratio = 1.0;
//! Default RT boundary sigma.
const float default_rt_sigma = 2;

/*! @brief Default minimum number of samples in score regions.
 *
 * Calculated from the default RT width and default RT sigma.
 */
const double root2pi = sqrt(2.0 * M_PI);
// The name of the program
const std::string program_name = "HiTIME";
// Number of spectra to keep in cache when reading from file.
// If the number is too small we may end up reading the same spectrum from disk multiple times
// and duplicating it in RAM. This may also limit multithreaded scalability.
// If we make it too big we might keep spectra in memory longer than needed and thus may
// waste some space.
const int default_input_spectrum_cache_size = 50;

#endif
