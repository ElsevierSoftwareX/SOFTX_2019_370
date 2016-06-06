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
   double_vect score;
   for (int n = low_spectra; n < high_spectra; n++)
   {
       score = score_spectra(input_map, n, half_window, opts);

       MSSpectrum<> input_spectrum = input_map.getSpectrum(n);
       MSSpectrum<> output_spectrum = MSSpectrum<Peak1D>(input_spectrum);

       // Copy the computed score into the output intensity
       for (int index = 0; index < input_spectrum.size(); ++index)
       {
           output_spectrum[index].setIntensity(score[index]);
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


/*! Calculate correlation scores for each MZ point in a central spectrum of
 * a data window.
 *
 * @param map collection of spectra in the input.
 * @param centre_idx Index of the spectrum to score.
 * @param half_window The number of spectra each side of the central spectrum
 * to include.
 * @param opts Options object.
 *
 * @return Vector of five vectors (min score, correlAB, correlA0, correlB0,
 * correl1r) giving score for each MZ in the central spectrum.
 */

double_vect
score_spectra(MSExperiment<> &map, int centre_idx, int half_window, Options opts)
{
    // Calculate constant values
    double rt_width_opt = opts.rt_width;
    double mz_width_opt = opts.mz_width;
    double mz_sigma_opt = opts.mz_sigma;
    double mz_delta_opt = opts.mz_delta;
    double min_sample_opt = opts.min_sample;
    double intensity_ratio_opt = opts.intensity_ratio;
    double rt_sigma = rt_width_opt / std_dev_in_fwhm;
    double mz_ppm_sigma = mz_width_opt / (std_dev_in_fwhm * 1e6);
    Size rt_len = map.getNrSpectra();
    // XXX why rename centre_idx to mid_win?
    int mid_win = centre_idx;
    double lo_tol = 1.0 - mz_sigma_opt * mz_ppm_sigma;
    double hi_tol = 1.0 + mz_sigma_opt * mz_ppm_sigma;
    int rt_offset = mid_win - half_window;

    double_vect score;

    //std::cerr << "Here 1\n";

    // XXX mz_mu_vect seems like a bad name
    MSSpectrum<> mz_mu_vect = map.getSpectrum(mid_win);

    // Low peak tolerances
    double_vect points_lo_lo;
    double_vect points_lo_hi;

    double_2d data_lo;
    std::vector<int> len_lo;

    //std::cerr << "Here 2\n";

    // Calculate tolerances for the lo and hi peak for each central MZ
    MSSpectrum<>::Iterator it;
    for (it = mz_mu_vect.begin(); it != mz_mu_vect.end(); ++it)
    {
	double this_mz = it->getMZ();
        points_lo_lo.push_back(this_mz * lo_tol);
        points_lo_hi.push_back(this_mz * hi_tol);

        double_vect data;
        data_lo.push_back(data);
        len_lo.push_back(0);
    }

    //std::cerr << "Here 3\n";

    // Iterate over the spectra in the window
    // Compute all the local data for each mz in the central spectrum
    for (int rowi = mid_win - half_window; rowi <= mid_win + half_window; ++rowi)
    {
	MSSpectrum<> rowi_spectrum;
        if (rowi >= 0 && rowi < rt_len)
	{
       	    rowi_spectrum = map.getSpectrum(rowi);
	}
        //XXX what do we do if this is false?

        // Iterate over the points in the central spectrum
        for (size_t mzi = 0; mzi < mz_mu_vect.size(); ++mzi)
	{
            // Get the tolerances and value for this point
            double lo_tol_lo = points_lo_lo[mzi];
            double lo_tol_hi = points_lo_hi[mzi];

            // Check if spectrum within bounds of the file...
            if (rowi >= 0 && rowi < rt_len)
	    {
                // Select points within tolerance for current spectrum
		Size lo_index = rowi_spectrum.findNearest(lo_tol_lo);
		Size hi_index = rowi_spectrum.findNearest(lo_tol_hi);

                // Check if points found...
		// XXX should check if the value found at index is near to our target mz
                if (lo_index <= hi_index)
		{
                    // Calculate guassian value for each found MZ
                    for (Size index = lo_index; index <= hi_index; ++index)
		    {
			Peak1D peak = rowi_spectrum[index];
		        double intensity = peak.getIntensity();
                        data_lo[mzi].push_back(intensity);
                    }
                    len_lo[mzi] += hi_index - lo_index + 1;

                // ...if not, use dummy data
                }
		else
		{
                    data_lo[mzi].push_back(0);
		    // XXX I think this should be:
		    len_lo[mzi] += 1;
                }

            // ...if outside use dummy data for this spectrum
            }
	    else
	    {
                data_lo[mzi].push_back(0);
		// XXX I think this should be:
		len_lo[mzi] += 1;
            }
        }
    }

    //std::cerr << "Here 4\n";
    
    // Check if the center intensity is a local maxima
    for (size_t mzi = 0; mzi < mz_mu_vect.size(); ++mzi)
    {
        //std::cerr << "Here 5\n";
	double center_intensity = mz_mu_vect[mzi].getIntensity();
        double result = center_intensity;
        for (int data_index = 0; data_index < len_lo[mzi]; ++data_index)
        {
           //std::cerr << "Here 6\n";
           if (center_intensity < data_lo[mzi][data_index])
           {
              //std::cerr << "Here 7\n";
              result = 0.0;
              break;
           }
        }
        //std::cerr << "Here 8\n";
        score.push_back(result);
    }

    //XXX fixme
    return score;
}

/*! Scores are output in CSV format with the following fields: retention time,
 * mz, intensity, minimum score, correlAB, correlA0, correlB0, correl1r.
 *
 * @param scores 2D vector of scores returned by score_spectra.
 * @param raw_data Spectrum pointer to the raw central vector.
 * @param out_stream Stream to write output too.
 * @param opts User defined Options object.
 */
void write_scores(double_2d scores, MSSpectrum<> raw_data, std::ofstream& out_stream)
{
    // Get central spectrum retention time
    double rt = raw_data.getRT();

    // Write output
    for (size_t idx = 0; idx < raw_data.size(); ++idx)
    {
        double mz  = raw_data[idx].getMZ();
        double amp = raw_data[idx].getIntensity();
        double ms  = scores[0][idx];
        double AB  = scores[1][idx];
        double A0  = scores[2][idx];
        double B0  = scores[3][idx];
        double r1  = scores[4][idx];

        out_stream << rt << ", " << mz << ", " << amp << ", "
                   << ms << ", " << AB << ", " << A0 << ", "
                   << B0 << ", " << r1 << std::endl;
    }
}
