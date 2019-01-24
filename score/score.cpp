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


Scorer::Scorer(bool debug, bool list_max, double intensity_ratio, double rt_width, 
               double mz_width, double mz_delta,
               double confidence,
               int num_threads, int input_spectrum_cache_size, string in_file, string out_file)
   : debug(debug)
   , list_max(list_max)
   , intensity_ratio(intensity_ratio)
   , rt_width(rt_width)
   , mz_width(mz_width)
   , mz_delta(mz_delta)
   , num_threads(num_threads)
   , confidence(confidence)
   , in_file(in_file)
   , out_file(out_file)
   , input_spectrum_cache(input_spectrum_cache_size)
   , spectrum_writer(out_file)
   , current_spectrum_id{0}
   , next_output_spectrum_id{0}
{

   IndexedMzMLFileLoader mzml;
   rt_sigma = default_rt_sigma;
   mz_sigma = default_mz_sigma;

   mzml.load(in_file, input_map);

   half_window = ceil(rt_sigma * rt_width / std_dev_in_fwhm);
   num_spectra = input_map.getNrSpectra();
   local_rows = (2 * half_window) + 1;

   Size min_sample = half_window;

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

       if (list_max)
       {
           score = local_max_spectra(this_spectrum_id);
       }
       else
       {
           score = score_spectra(this_spectrum_id);
       }
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

void Scorer::collect_local_rows(int rt_offset, double_2d &mz_vals, double_2d &amp_vals)
{
    PeakSpectrumPtr rowi_spectrum;
    // Iterate over the spectra in the window
    for (Size rowi = 0; rowi < mz_vals.size(); ++rowi)
    {
        // window can go outside start and end of scans, so check bounds
        if (rowi + rt_offset >= 0 && rowi + rt_offset < num_spectra)
        {
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

void Scorer::collect_window_data(double scale,
               double_vect & rt_shape,
               double centre, double sigma,
               double_2d & mz_vals, double_2d & amp_vals,
               double lower_bound_mz, double upper_bound_mz,
               double_vect & data_out, double_vect & shape_out)
{
    // Iterate over the spectra in the window
    for (Size rowi = 0; rowi < mz_vals.size() && rowi < rt_shape.size(); ++rowi)
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

        // Calculate Gaussian value for each found MZ
        for (Size index = lower_index; index <= upper_index && index < mz_vals[rowi].size(); ++index)
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

double correlation(double_vect &data, double_vect &shape)
{
    // Zero correlation if not enough data in either region
    if (data.size() < 2 or shape.size() < 2) return 0.0;

    // local copies
    double_vect X (data);
    double_vect Y (shape);

    double EX = std::accumulate(X.begin(), X.end(), 0.0) / X.size();
    double EY = std::accumulate(Y.begin(), Y.end(), 0.0) / Y.size();

    // Centre data in regions relative to combined means
    for (auto& val : X) val -= EX;
    for (auto& val : Y) val -= EY;

    // COV and VAR
    double COV = inner_product(X.begin(), X.end(), Y.begin(), 0.0);
    double VARX = inner_product(X.begin(), X.end(), X.begin(), 0.0);
    double VARY = inner_product(Y.begin(), Y.end(), Y.begin(), 0.0);

    double correl = COV / std::sqrt(VARX * VARY);
    if (std::isnan(correl) or std::isinf(correl)) correl = 0.0;

    return correl;
}

double combined_correlation(double_vect &data_nat, double_vect &data_iso,
                double_vect &shape_nat, double_vect &shape_iso)
{
    // local copies
    double_vect Xa (data_nat);
    double_vect Xb (data_iso);
    double_vect Ya (shape_nat);
    double_vect Yb (shape_iso);

    double EXa = std::accumulate(Xa.begin(), Xa.end(), 0.0) / Xa.size();
    double EXb = std::accumulate(Xb.begin(), Xb.end(), 0.0) / Xb.size();
    double EYa = std::accumulate(Ya.begin(), Ya.end(), 0.0) / Ya.size();
    double EYb = std::accumulate(Yb.begin(), Yb.end(), 0.0) / Yb.size();

    // Combined mean is mean of means
    double EXab = 0.5 * (EXa + EXb);  // E(Xab) =  E(E(Xa), E(Xb))
    double EYab = 0.5 * (EYa + EYb);  // E(Yab) =  E(E(Ya), E(Yb))

    // Centre data in regions relative to combined means
    for (auto& val : Xa) val -= EXab;  // (Xa - E(Xab)) --> C(Xa) 
    for (auto& val : Xb) val -= EXab;  // (Xb - E(Xab))
    for (auto& val : Ya) val -= EYab;  // (Ya - E(Yab))
    for (auto& val : Yb) val -= EYab;  // (Yb - E(Yab))

    // region expected values
    // E((Xa - E(Xab))(Ya - E(Yab)))
    // E((Xa - E(Xab))^2)
    // E((Ya - E(Yab))^2)
    double ECXaCYa = inner_product(Xa.begin(), Xa.end(), Ya.begin(), 0.0);
    double ECXa2 = inner_product(Xa.begin(), Xa.end(), Xa.begin(), 0.0);
    double ECXb2 = inner_product(Xb.begin(), Xb.end(), Xb.begin(), 0.0);

    // E((Xb - E(Xab))(Yb - E(Yab)))
    // E((Xb - E(Xab))^2)
    // E((Yb - E(Yab))^2)
    double ECXbCYb = inner_product(Xb.begin(), Xb.end(), Yb.begin(), 0.0);
    double ECYa2 = inner_product(Ya.begin(), Ya.end(), Ya.begin(), 0.0);
    double ECYb2 = inner_product(Yb.begin(), Yb.end(), Yb.begin(), 0.0);

    // Combined COV, VAR as mean region COV and VAR
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
double mengZ(double rhoXY, double rhoXZ, double rhoYZ, Size samples, double confidence)
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
        if (confidence > 0.0) {
            // Only use score if lower confidence greater than zero
            double lCI = (std::atanh(rhoXY) - std::atanh(rhoXZ)) - confidence *
                    std::sqrt((2.0 * (1.0 - rhoYZ) * h) / (samples - 3.0));
            if (std::isnan(lCI) or std::isinf(lCI) or lCI < 0.0) lCI = 0.0;
            if (lCI > 0.0) return z;
        }

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
    double lower_tol = 1.0 - mz_sigma * mz_ppm_sigma;
    double upper_tol = 1.0 + mz_sigma * mz_ppm_sigma;
    int rt_offset = centre_idx - half_window;

    // XXX This should be calculated once, then used as look-up
    // Spacing should also be based on scan intervals
    // (curently assumed fixed spacing)
    double_vect rt_shape (local_rows);

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
    collect_local_rows(rt_offset, mz_vals, amp_vals);

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

        collect_window_data(1.0,  rt_shape,
                        centre, sigma, mz_vals, amp_vals,
                        lower_bound_nat, upper_bound_nat, data_nat, shape_nat);
        collect_window_data(intensity_ratio, rt_shape,
                        centre_iso, sigma_iso, mz_vals, amp_vals,
                        lower_bound_iso, upper_bound_iso, data_iso, shape_iso);

        // Zero score if not enough data in either region
        if (data_nat.size() < min_sample or data_iso.size() < min_sample)
        {
            min_score_vect.push_back(0.0);
            continue;
        }

        /*
         * Competing models
         * Target model with desired isotope ion ratio
         * Model 1, higher ratio
         * Model 2, lower ratio
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

        // Only contrast if natural ion correlates to model
        // User lower confidence interval at given confidence
        if (confidence > 0.0) {
            double z1 = correlation(data_nat, shape_nat);
            z1 = std::atanh(z1) - confidence/std::sqrt(data_nat.size() - 3.0);
            if (std::isnan(z1) or std::isinf(z1) or z1 <= 0.0)
            {
                min_score_vect.push_back(0.0);
                continue;
            }
            // Only contrast if isotope ion correlates to model
            z1 = correlation(data_iso, shape_iso);
            z1 = std::atanh(z1) - confidence/std::sqrt(data_nat.size() - 3.0);
            if (std::isnan(z1) or std::isinf(z1) or z1 <= 0.0)
            {
                min_score_vect.push_back(0.0);
                continue;
            }
        }

        /* Alternate models */
        // Twin ion with different ratios
        double correl_XabYab = combined_correlation(data_nat, data_iso, shape_nat, shape_iso);
        double_vect iso_lower (shape_iso);
        for (auto& val : iso_lower) val *= 0.001;
        double correl_XabYa_ = combined_correlation(data_nat, data_iso, shape_nat, iso_lower);
        // inter shape correlation
        double correl_YabYa_ = combined_correlation(shape_nat, shape_iso, shape_nat, iso_lower);
 
        double_vect nat_lower (shape_nat);
        for (auto& val : nat_lower) val *= 0.001;
        double correl_XabY_b = combined_correlation(data_nat, data_iso, nat_lower, shape_iso);
        // inter shape correlation
        double correl_YabY_b = combined_correlation(shape_nat, shape_iso, nat_lower, shape_iso);

        // Calculate z scores
        nAB = data_nat.size() + data_iso.size();
        double zABA0 = mengZ(correl_XabYab, correl_XabYa_, correl_YabYa_, nAB, confidence);
        double zAB0B = mengZ(correl_XabYab, correl_XabY_b, correl_YabY_b, nAB, confidence);
    
        // Find the minimum scores, bounded at zero
        double min_score = std::max({0.0, std::min({zABA0, zAB0B})});
    
        min_score_vect.push_back(min_score);
    } 
    return min_score_vect;
}


bool Scorer::local_max_data(double centre_amp,
               double_2d & mz_vals, double_2d & amp_vals,
               double lower_bound_mz, double upper_bound_mz)
{
    // Iterate over the spectra in the window
    for (Size rowi = 0; rowi < mz_vals.size(); ++rowi)
    {
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

        // Calculate Gaussian value for each found MZ
        for (Size index = lower_index; index <= upper_index && index < mz_vals[rowi].size(); ++index)
        {
            double mz = mz_vals[rowi][index];
            double intensity = amp_vals[rowi][index];

            // just in case
            if (mz < lower_bound_mz || mz > upper_bound_mz) continue;

            if (intensity > centre_amp)
                return false;

        }
    }
    return true;
}


double_vect Scorer::local_max_spectra(int centre_idx)
{
    // Calculate constant values
    double local_rt_sigma = rt_width / std_dev_in_fwhm;
    int rt_offset = centre_idx - half_window;

    PeakSpectrumPtr centre_row_points = get_spectrum(centre_idx);

    // Length of all vectors (= # windows)
    Size mz_windows = centre_row_points->size();

    double_vect min_score_vect;
    min_score_vect.reserve(mz_windows);

    // NOTE: much faster to collect all local row data 1st
    double_2d mz_vals (local_rows);
    double_2d amp_vals (local_rows);

    // collect all row data
    collect_local_rows(rt_offset, mz_vals, amp_vals);

    PeakSpectrum::Iterator it;
    for (it = centre_row_points->begin(); it != centre_row_points->end(); ++it)
    {
        double centre_mz = it->getMZ();
        double centre_amp = it->getIntensity();

        if (local_max_data(centre_amp, mz_vals, amp_vals,
                           centre_mz - mz_width, centre_mz + mz_width))
        {
            min_score_vect.push_back(centre_amp);
        }
        else
        {
            min_score_vect.push_back(0.0);
        }
    } 
    return min_score_vect;
}
