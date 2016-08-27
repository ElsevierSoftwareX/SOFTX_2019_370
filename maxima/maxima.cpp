#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <mutex>
#include <iostream>
#include "vector.h"
#include "options.h"
#include "constants.h"
#include "maxima.h"

using namespace OpenMS;
using namespace std;

mutex output_spectrum_lock;

// Worker function for scoring a block of spectra. This is executed by a thread. 
void score_worker(MSExperiment<> &input_map, MSExperiment<> &output_map, int half_window, Options opts, int low_spectra, int high_spectra)
{
   for (int n = low_spectra; n < high_spectra; n++)
   {
       double_vect score = score_spectra(input_map, n, half_window, opts);
       MSSpectrum<> input_spectrum = input_map.getSpectrum(n);
       MSSpectrum<> output_spectrum = MSSpectrum<Peak1D>(input_spectrum);
       double this_rt = input_spectrum.getRT();
       double this_mz;
       double this_score;

       // Copy the computed score into the output intensity
       for (int index = 0; index < input_spectrum.size(); ++index)
       {
           output_spectrum[index].setIntensity(score[index]);
           this_score = output_spectrum[index].getIntensity();
           if (this_score > 0)
           {
               this_mz = output_spectrum[index].getMZ();
               //cout << this_rt << " " << this_mz << " " << this_score << endl;
               cout << fixed << this_mz << " " << fixed << this_rt << " " << fixed << this_score << endl;
           }
       }
       // Writing to the output_map must be synchronised between threads, with only one thread
       // allowed to write at a time.
       unique_lock<mutex> lck {output_spectrum_lock};
       output_map.addSpectrum(output_spectrum);
       output_spectrum_lock.unlock();
       // Output the results
       // write_scores(score, input_spectrum, outfile);
   }
}

/* Find index of the least mass >= low_mass, and the greatest mass <= high_mass */

void get_bounds(MSSpectrum<> &spectrum, double low_mass, double high_mass, Size &lo_index, Size &hi_index)
{
    double this_mass;
    lo_index = -1;
    hi_index = -1;
    Size probe;
    Size max_probe = spectrum.size() - 1;

    probe = spectrum.findNearest(low_mass);

    if (probe >= 0 && probe <= max_probe)
    {
       this_mass = spectrum[probe].getMZ();
       while(probe <= max_probe && this_mass < low_mass)
       {
          probe += 1;
          this_mass = spectrum[probe].getMZ();
       }

       if (this_mass >= low_mass && this_mass <= high_mass)
       {
          lo_index = probe;
       }

    }

    probe = spectrum.findNearest(high_mass);

    if (probe >= 0 && probe <= max_probe)
    {
       this_mass = spectrum[probe].getMZ();
       while(probe > 0 && this_mass > high_mass)
       {
          probe -= 1;
          this_mass = spectrum[probe].getMZ();
       }

       if (this_mass >= low_mass && this_mass <= high_mass)
       {
          hi_index = probe;
       }
    }
}



/*! Compute the local maxima scores for the input spectrum.
 *
 * @param map collection of spectra in the input.
 * @param centre_idx Index of the spectrum to score.
 * @param half_window The number of spectra each side of the central spectrum
 * to include.
 * @param opts Options object.
 *
 * @return Vector of scores for the central spectrum.
 */

double_vect
score_spectra(MSExperiment<> &map, int mid_win, int half_window, Options opts)
{
    double mz_width_opt = opts.mz_width;
    double mz_sigma_opt = opts.mz_sigma;
    double mz_ppm_sigma = mz_width_opt / (std_dev_in_fwhm * 1e6);
    Size rt_len = map.getNrSpectra();
    double lo_tol = 1.0 - mz_sigma_opt * mz_ppm_sigma;
    double hi_tol = 1.0 + mz_sigma_opt * mz_ppm_sigma;

    MSSpectrum<> mz_mu_vect = map.getSpectrum(mid_win);

    // Low peak tolerances
    double_vect points_lo;
    double_vect points_hi;
    double_2d local_intensities;

    // Calculate tolerances for the lo and hi peak for each central MZ
    MSSpectrum<>::Iterator it;
    for (it = mz_mu_vect.begin(); it != mz_mu_vect.end(); ++it)
    {
	double this_mz = it->getMZ();
        points_lo.push_back(this_mz * lo_tol);
        points_hi.push_back(this_mz * hi_tol);

        // XXX does this share the vector, or create a new one each time?
        double_vect data;
        local_intensities.push_back(data);
    }

    // Iterate over the spectra in the window
    // Compute all the local data for each mz in the central spectrum
    for (int rowi = mid_win - half_window; rowi <= mid_win + half_window; ++rowi)
    {
       if (rowi >= 0 && rowi < rt_len)
       {
          MSSpectrum<> rowi_spectrum = map.getSpectrum(rowi);
          // Iterate over the points in the central spectrum
          for (size_t mzi = 0; mzi < mz_mu_vect.size(); ++mzi)
	  {
             // Get the tolerances and value for this point
             double lo_tol_lo = points_lo[mzi];
             double lo_tol_hi = points_hi[mzi];

             // Select points within tolerance for current spectrum

             Size lo_index;
             Size hi_index;
             get_bounds(rowi_spectrum, lo_tol_lo, lo_tol_hi, lo_index, hi_index);
             // Check if points found...
             // XXX should check if the value found at index is near to our target mz
             if (lo_index != -1 && hi_index != -1 && lo_index <= hi_index)
             {
                // Collect neighbouring local intensities
                for (Size index = lo_index; index <= hi_index; ++index)
	        {
                   double intensity = rowi_spectrum[index].getIntensity();
                   local_intensities[mzi].push_back(intensity);
                }
             }
          }
       }
    }

    double_vect scores;

    // Check if the center intensity is a local maxima
    for (size_t mzi = 0; mzi < mz_mu_vect.size(); ++mzi)
    {
	double center_intensity = mz_mu_vect[mzi].getIntensity();
        double result = center_intensity;

        for (size_t data_index = 0; data_index < local_intensities[mzi].size(); ++data_index)
        {
           if (center_intensity < local_intensities[mzi][data_index])
           {
              result = 0.0;
              break;
           }
        }
        scores.push_back(result);
    }

    return scores;
}
