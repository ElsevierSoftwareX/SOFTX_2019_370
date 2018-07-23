#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <mutex>
#include <iostream>
#include <map>
#include <thread>
#include "vector.h"
#include "options.h"
#include "constants.h"
#include "lru_cache.h"
#include "score.h"

// input cache size
// #define CACHE_SIZE 30
// #define CACHE_SIZE 1000 

using namespace OpenMS;
using namespace std;

mutex output_spectrum_lock;
mutex next_spectrum_lock;
mutex input_spectrum_lock;


Scorer::Scorer(bool debug, double intensity_ratio, double rt_width, double rt_sigma, double ppm,
               double mz_width, double mz_sigma, double mz_delta, double min_sample, int num_threads,
               int input_spectrum_cache_size, string in_file, string out_file)
   : debug(debug)
   , intensity_ratio(intensity_ratio)
   , rt_width(rt_width)
   , rt_sigma(rt_sigma)
   , ppm(ppm)
   , mz_width(mz_width)
   , mz_sigma(mz_sigma)
   , mz_delta(mz_delta)
   , min_sample(min_sample)
   , num_threads(num_threads)
   , in_file(in_file)
   , out_file(out_file)
   , input_spectrum_cache(input_spectrum_cache_size)
   , spectrum_writer(out_file)
   , current_spectrum_id{0}
   , next_output_spectrum_id{0}
{

   IndexedMzMLFileLoader mzml;

   mzml.load(in_file, input_map);

   half_window = ceil(rt_sigma * rt_width / std_dev_in_fwhm);
   num_spectra = input_map.getNrSpectra();

   vector<thread> threads(num_threads);

   //cout << "Num threads: " << num_threads << endl;
   //cout << "Num spectra: " << num_spectra << endl;

   for (int thread_count = 0; thread_count < num_threads; thread_count++)
   {
      threads[thread_count] = thread(&Scorer::score_worker, this, thread_count);
    
   }

   for (int thread_count = 0; thread_count < num_threads; thread_count++)
   {
       threads[thread_count].join();
   }
}

/*
PeakSpectrumPtr Scorer::get_spectrum(int spectrum_id)
{
   PeakSpectrumPtr spectrum_ptr;

   input_spectrum_lock.lock();
   spectrum_ptr = make_shared<PeakSpectrum>(input_map.getSpectrum(spectrum_id));
   input_spectrum_lock.unlock();
   return spectrum_ptr;
}
*/
PeakSpectrumPtr Scorer::get_spectrum(int spectrum_id)
{
   PeakSpectrumPtr spectrum_ptr;

   input_spectrum_lock.lock();
   if (input_spectrum_cache.exists(spectrum_id))
   {
      spectrum_ptr = input_spectrum_cache.get(spectrum_id);
   }
   else
   {
      spectrum_ptr = make_shared<PeakSpectrum>(input_map.getSpectrum(spectrum_id));
      input_spectrum_cache.put(spectrum_id, spectrum_ptr);
   }
   input_spectrum_lock.unlock();
   return spectrum_ptr;
}

void Scorer::put_spectrum(int spectrum_id, PeakSpectrum spectrum)
{
   output_spectrum_lock.lock();

   if (spectrum_id == next_output_spectrum_id)
   {
      // this is the next spectrum to output
      spectrum_writer.consumeSpectrum(spectrum);
      next_output_spectrum_id++;

      // try to output more spectra
      while(output_spectrum_queue.size() > 0)
      {
         IndexSpectrum index_spectrum = output_spectrum_queue.top();
         if (index_spectrum.first == next_output_spectrum_id)
         {
            spectrum_writer.consumeSpectrum(index_spectrum.second);
            output_spectrum_queue.pop();
            next_output_spectrum_id++;
         }
         else
         {
            break;
         }
      }
   }
   else
   {
      // push this spectrum into the queue to write out later
      output_spectrum_queue.push(pair<int, PeakSpectrum>(spectrum_id, spectrum));
   } 

   output_spectrum_lock.unlock();
}

int Scorer::get_next_spectrum_todo(void)
{
   int this_spectrum;
   next_spectrum_lock.lock();
   this_spectrum = current_spectrum_id;
   current_spectrum_id++;
   next_spectrum_lock.unlock();
   return this_spectrum; 
}

void Scorer::score_worker(int thread_count)
{
   double_vect score;
   int this_spectrum_id;

   this_spectrum_id = get_next_spectrum_todo(); 

   while (this_spectrum_id < num_spectra)
   {
       //if ((this_spectrum_id % 100) == 0) {
       //    cout << "Thread: " << thread_count << " Spectrum: " << this_spectrum_id << endl;
       //}

       score = score_spectra(this_spectrum_id);
       PeakSpectrumPtr input_spectrum = get_spectrum(this_spectrum_id);
       
       PeakSpectrum output_spectrum = MSSpectrum<Peak1D>(*input_spectrum);
       for (int index = 0; index < input_spectrum->size(); index++)
       {
           output_spectrum[index].setIntensity(score[index]);
       }

       put_spectrum(this_spectrum_id, output_spectrum);
       
       this_spectrum_id = get_next_spectrum_todo(); 
   }

}

void Scorer::collect_local_rows(Size rt_offset, Size rt_len,
                double_2d &mz_vals, double_2d &amp_vals)
{
    // Iterate over the spectra in the window
    for (Size rowi = 0; rowi < mz_vals.size(); ++rowi)
    {
        // window can go outside start and end of scans, so check bounds
        if (rowi + rt_offset >= 0 && rowi + rt_offset < rt_len)
        {
            PeakSpectrumPtr rowi_spectrum;
            rowi_spectrum = get_spectrum(rowi + rt_offset);
            Size elements = rowi_spectrum->size();

            mz_vals[rowi].resize(elements);
            amp_vals[rowi].resize(elements);

            // This shouldn't happen, so want to but want know if it does
            // XXX maybe we should provide an option to avoid this for performance reasons?
            if (!rowi_spectrum->isSorted()) throw std::runtime_error ("Spectrum not sorted");

            for (Size index = 0; index < elements; ++index)
            {
                Peak1D peak = (*rowi_spectrum)[index];
                double mz = peak.getMZ();
                double intensity = peak.getIntensity();

                mz_vals[rowi][index] = mz;
                amp_vals[rowi][index] = intensity;
            }
        }
        else
        {
            mz_vals[rowi].clear();
            amp_vals[rowi].clear();
        }
    }
}

void Scorer::collect_window_data(Size rt_offset, Size rt_len, double scale,
               double_vect & rt_shape,
               double centre, double sigma,
               double_2d & mz_vals, double_2d & amp_vals,
               double lower_bound_mz, double upper_bound_mz,
               double_vect & data_out, double_vect & shape_out)
{
        // Iterate over the spectra in the window
        for (Size rowi = 0; rowi < mz_vals.size(); ++rowi)
        {
            // window can go outside start and end of scans, so check bounds
            if (rowi + rt_offset >= 0 && rowi + rt_offset < rt_len)
            {
                double rt_shape_i = rt_shape[rowi] * scale;

                // Select points within tolerance for current spectrum
                // Want index of bounds
                // Need to convert iterator to index
                Size lower_index = Size(std::lower_bound(mz_vals[rowi].begin(),
                                            mz_vals[rowi].end(),
                                            lower_bound_mz)
                                        -
                                        mz_vals[rowi].begin());
                if (lower_index < 0) lower_index = 0;
                Size upper_index = Size(std::lower_bound(mz_vals[rowi].begin(),
                                            mz_vals[rowi].end(),
                                            upper_bound_mz)
                                        -
                                        mz_vals[rowi].begin());
                // this test should also stop lookup in empty rows
                if (upper_index >= mz_vals[rowi].size()) upper_index = mz_vals[rowi].size() - 1;

                // Calculate Gaussian value for each found MZ
                for (Size index = lower_index; index <= upper_index; ++index)
                {
                    double mz = mz_vals[rowi][index];
                    double intensity = amp_vals[rowi][index];

                    // just in case
                    if (mz < lower_bound_mz || mz > upper_bound_mz) continue;

                    // calc mz fit
                    mz = (mz - centre) / sigma;
                    mz = -0.5 * mz * mz;
                    double fit = exp(mz) / (sigma * root2pi);

                    data_out.push_back(intensity);
                    shape_out.push_back(fit * rt_shape_i);
                }
            }
        }
}

double combined_correlation(double_vect &data_nat, double_vect &data_iso,
                double_vect &shape_nat, double_vect &shape_iso)
{
    // Region means
    // E(Xa) 
    // E(Xb)
    // E(Ya)
    // E(Yb)
    double EXa = mean_vector(data_nat);
    double EXb = mean_vector(data_iso);
    double EYa = mean_vector(shape_nat);
    double EYb = mean_vector(shape_iso);

    // Combined mean is mean of means
    // E(Xab) =  E(E(Xa), E(Xb))
    // E(Yab) =  E(E(Ya), E(Yb))
    double EXab = 0.5 * (EXa + EXb);
    double EYab = 0.5 * (EYa + EYb);

    // Centre data in regions relative to combined means
    double_vect CXa (data_nat);    // (Xa - E(Xab)) --> C(Xa)
    double_vect CXb (data_iso);    // (Xb - E(Xab))
    double_vect CYa (shape_nat);    // (Ya - E(Yab))
    double_vect CYb (shape_iso);    // (Yb - E(Yab))

    // Centre data in regions relative to combined means
    for (auto& val : CXa) val -= EXab;
    for (auto& val : CXb) val -= EXab;
    for (auto& val : CYa) val -= EYab;
    for (auto& val : CYb) val -= EYab;

    // region expected values
    // E((Xa - E(Xab))(Ya - E(Yab)))
    // E((Xa - E(Xab))^2)
    // E((Ya - E(Yab))^2)

    // E((Xb - E(Xab))(Yb - E(Yab)))
    // E((Xb - E(Xab))^2)
    // E((Yb - E(Yab))^2)

    double ECXaCYa = inner_product(CXa.begin(), CXa.end(), CYa.begin(), 0.0);
    double ECXbCYb = inner_product(CXb.begin(), CXb.end(), CYb.begin(), 0.0);
    double ECXa2 = inner_product(CXa.begin(), CXa.end(), CXa.begin(), 0.0);
    double ECXb2 = inner_product(CXb.begin(), CXb.end(), CXb.begin(), 0.0);
    double ECYa2 = inner_product(CYa.begin(), CYa.end(), CYa.begin(), 0.0);
    double ECYb2 = inner_product(CYb.begin(), CYb.end(), CYb.begin(), 0.0);

    double cov_Xab = 0.5 * (ECXaCYa + ECXbCYb);
    double var_Xab = 0.5 * (ECXa2 + ECXb2);
    double var_Yab = 0.5 * (ECYa2 + ECYb2);

    double correl = cov_Xab / std::sqrt(var_Xab * var_Yab);
    if (std::isnan(correl) or std::isinf(correl)) correl = 0.0;

    return correl;
}

/*
 * Calculate Meng's Z-score
 * Meng, Rubin, & Rosenthal (1992),
 * Comparing Correlated Correlation Coefficients,
 * Psychological Bulletin 111(1), 172-175.
 */
double Scorer::mengZ(double rhoXY, double rhoXZ, double rhoYZ, Size samples)
{
        // Calculate rm values between correlations
        double rm2 = 0.5 * (rhoXY * rhoXY + rhoXZ * rhoXZ);
    
        // Calculate f values between correlation and rm
        double f = (1.0 - rhoYZ) / (2.0 * (1.0 - rm2));
        if (std::isinf(f) or f > 1.0) f = 1.0;
    
        // Calculate h values between f and rm
        double h = (1.0 - f * rm2) / (1.0 - rm2);
    
        // Calculate z scores
        double z = (std::atanh(rhoXY) - std::atanh(rhoXZ)) *
               std::sqrt( (samples - 3.0) / (2.0 * (1.0 - rhoYZ) * h) );
        if (std::isnan(z) or std::isinf(z) or z < 0.0) z = 0.0;

        // only return Z score if lower confidence interval > 0
        // TODO: make CONFIDENCE a parameter
        double CONFIDENCE = 1.96;  // 5% error

        return z;
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
 * @return Vector of score at each MZ in central spectrum.
 */

double_vect Scorer::score_spectra(int centre_idx)
{
    // Calculate constant values
    double local_rt_sigma = rt_width / std_dev_in_fwhm;
    double mz_ppm_sigma = mz_width / (std_dev_in_fwhm * 1e6);
    Size rt_len = input_map.getNrSpectra();
    Size local_rows = (2 * half_window) + 1;
    double lower_tol = 1.0 - mz_sigma * mz_ppm_sigma;
    double upper_tol = 1.0 + mz_sigma * mz_ppm_sigma;
    Size rt_offset = centre_idx - half_window;

    // XXX This should be calculated once, then used as look-up
    // Spacing should also be based on scan intervals
    // (curently assumed fixed spacing)
    double_vect rt_shape (local_rows);
//    rt_shape.resize(local_rows);

    // Calculate Gaussian shape in the RT direction
    for (Size i = 0; i < local_rows; ++i)
    {
        double pt = (i - half_window) / local_rt_sigma;
        pt = -0.5 * pt * pt;
        double fit = exp(pt) / (local_rt_sigma * root2pi);
        rt_shape[i] = fit;
    }

    PeakSpectrumPtr centre_row_points = get_spectrum(centre_idx);

    // Length of all vectors (= # windows)
    Size mz_windows = centre_row_points->size();

    double_vect min_score_vect;
    min_score_vect.reserve(mz_windows);

    // Low (natural ion) peak tolerances
    double lower_bound_nat = 0.0;
    double upper_bound_nat = 0.0;
    // High (isotope ion) peak tolerances
    double lower_bound_iso = 0.0;
    double upper_bound_iso = 0.0;

    // total data per rt,mz centre
    double nAB = 0.0;

    double centre = 0.0;
    double sigma = 0.0;
    double centre_iso = 0.0;
    double sigma_iso = 0.0;

    // NOTE: much faster to collect all local row data 1st
    double_2d mz_vals (local_rows);
    double_2d amp_vals (local_rows);

    // collect all row data
    collect_local_rows(rt_offset, rt_len, mz_vals, amp_vals);

    // Calculate tolerances for the lo and hi peak for each central MZ
    double_vect data_nat;
    double_vect data_iso;
    double_vect shape_nat;
    double_vect shape_iso;

    PeakSpectrum::Iterator it;
    for (it = centre_row_points->begin(); it != centre_row_points->end(); ++it)
    {
        centre = it->getMZ();
        sigma = centre * mz_ppm_sigma;

        centre_iso = centre + mz_delta;
        sigma_iso = centre_iso * mz_ppm_sigma;

        lower_bound_nat = centre * lower_tol;
        upper_bound_nat = centre * upper_tol;
        lower_bound_iso = centre_iso * lower_tol;
        upper_bound_iso = centre_iso * upper_tol;

        // reset index back to start
        data_nat.clear();
        data_iso.clear();
        shape_nat.clear();
        shape_iso.clear();

        collect_window_data(rt_offset, rt_len, 1.0,  rt_shape,
                        centre, sigma, mz_vals, amp_vals,
                        lower_bound_nat, upper_bound_nat, data_nat, shape_nat);
        collect_window_data(rt_offset, rt_len, intensity_ratio, rt_shape,
                        centre_iso, sigma_iso, mz_vals, amp_vals,
                        lower_bound_iso, upper_bound_iso, data_iso, shape_iso);

        // Ignore regions with insufficient number of samples
        // If any region is ignored, set all to empty
        if (data_nat.size() < min_sample || data_iso.size() < min_sample)
        {
            data_nat.clear();
            shape_nat.clear();
            data_iso.clear();
            shape_iso.clear();
            nAB = 0;
        }
        else
        {
            nAB = data_nat.size() + data_iso.size();
        }

        /*
         * Competing models
         * Low ion window, High ion window
         * Empty low ion window, flat high ion window
         */

        /* Formulation
         * Correlation based on expectations in each region
         * Low ion region, a
         * High ion region, b
         * Correlation is
         * Covariance = E(E((Xa - E(Xab))(Ya - E(Yab))), E((Xb - E(Xab))(Yb -E(Yab))))
         * Data Variance = E(E((Xa - E(Xab))^2), E((Xb - E(Xab))^2))
         * Model Variance = E(E((Ya - E(Yab))^2), E((Yb - E(Yab))^2))
         */

        /* Alternate models */
        // low region and high region modelled as flat

        // Correlations
        double correl_XabYab = combined_correlation(data_nat, data_iso, shape_nat, shape_iso);

        double_vect shape_flat (shape_iso.size(), 0.0);
        double correl_XabYa_ = combined_correlation(data_nat, data_iso, shape_nat, shape_flat);

        // Yab, Ya_ covariance...
        double correl_YabYa_ = combined_correlation(shape_nat, shape_iso, shape_nat, shape_flat);

        shape_flat.clear();
        shape_flat.resize(shape_nat.size(), 0.0);
        double correl_XabY_b = combined_correlation(data_nat, data_iso, shape_flat, shape_iso);

        // Yab, Y_b covariance...
        double correl_YabY_b = combined_correlation(shape_nat, shape_iso, shape_flat, shape_iso);

        // Calculate z scores
        double zABA0 = mengZ(correl_XabYab, correl_XabYa_, correl_YabYa_, nAB);
        double zAB0B = mengZ(correl_XabYab, correl_XabY_b, correl_YabY_b, nAB);
    
        // Find the minimum scores, bounded at zero
        double min_score = std::max({0.0, std::min({zABA0, zAB0B})});
    
        min_score_vect.push_back(min_score);
    } 
    return min_score_vect;
}
