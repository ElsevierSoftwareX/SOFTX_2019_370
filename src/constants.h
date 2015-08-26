#include <math.h>

//! Convert standard deviation to FWHM
const float std_dev_in_fwhm = 2.355;
//! Default difference in mass of isotopes.
const float default_mz_delta        = 6.0201;
//! Default MZ tolerance in parts per million.
const float default_ppm             = 4.0;
//! Defualt MZ Full Width Half Maximum in PPM.
const float default_fwhm            = 150.0;
//! Default MZ boundary sigma.
const float default_mz_sigma        = 1.5;
//! Default ratio of peak intensities.
const float default_intensity_ratio = 1.0;
//! Default retention time FWHM in scans.
const float default_rt_width        = 17.0;
//! Default RT boundary sigma.
const float default_rt_sigma        = 1.5;
/*! @brief Default minimum number of samples in score regions.
 *
 * Calculated from the default RT width and default RT sigma.
 */
const float default_min_sample = default_rt_width * default_rt_sigma / std_dev_in_fwhm;
const double root2pi = sqrt(2.0 * M_PI);
