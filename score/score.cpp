#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <mutex>
#include <iostream>
#include "vector.h"
#include "options.h"
#include "constants.h"
#include "score.h"

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
 * @return Vector of score at each MZ in central spectrum.
 */

double_vect
score_spectra(MSExperiment<> &map, int centre_idx, int half_window, Options opts)
{
    // Calculate constant values
    double mz_delta_opt = opts.mz_delta;
    double min_sample_opt = opts.min_sample;
    double intensity_ratio_opt = opts.intensity_ratio;
    double rt_sigma = opts.rt_width / std_dev_in_fwhm;
    double mz_ppm_sigma = opts.mz_width / (std_dev_in_fwhm * 1e6);
    Size rt_len = map.getNrSpectra();
    double lower_tol = 1.0 - opts.mz_sigma * mz_ppm_sigma;
    double upper_tol = 1.0 + opts.mz_sigma * mz_ppm_sigma;
    int rt_offset = centre_idx - half_window;

    MSSpectrum<> centre_row_points = map.getSpectrum(centre_idx);

    // Length of all vectors (= # windows)
    size_t mz_windows = centre_row_points.size();

    // Low (natural ion) peak tolerances
    double_vect lower_bound_nat;
    double_vect upper_bound_nat;
    // High (isotope ion) peak tolerances
    double_vect lower_bound_iso;
    double_vect upper_bound_iso;

    // total data per rt,mz centre
    double_vect nAB;

    double_2d data_nat;
    double_2d data_iso;
    double_2d shape_nat;
    double_2d shape_iso;

    // Calculate tolerances for the lo and hi peak for each central MZ
    MSSpectrum<>::Iterator it;
    for (it = centre_row_points.begin(); it != centre_row_points.end(); ++it)
    {
        double this_mz = it->getMZ();
        lower_bound_nat.push_back(this_mz * lower_tol);
        upper_bound_nat.push_back(this_mz * upper_tol);
        lower_bound_iso.push_back((this_mz + mz_delta_opt) * lower_tol);
        upper_bound_iso.push_back((this_mz + mz_delta_opt) * upper_tol);

        double_vect data;
        data_nat.push_back(data);
        data_iso.push_back(data);
        shape_nat.push_back(data);
        shape_iso.push_back(data);
    }

    // XXX This should be calculated once, then used as look-up
    std::vector<double> rt_shape;

    // Calculate Gaussian shape in the RT direction
    for (int i = 0; i < (2 * half_window) + 1; ++i)
    {
        double pt = (i - half_window) / rt_sigma;
        pt = -0.5 * pt * pt;
        double fit = exp(pt) / (rt_sigma * root2pi);
        rt_shape.push_back(fit);
    }

    // Iterate over the spectra in the window
    for (int rowi = centre_idx - half_window; rowi <= centre_idx + half_window; ++rowi)
    {
        double rt_shape_nat = rt_shape[rowi - rt_offset];
        double rt_shape_iso = rt_shape_nat * intensity_ratio_opt;

        MSSpectrum<> rowi_spectrum;
        // window can go outside start and end of scans, so check bounds
        if (rowi >= 0 && rowi < rt_len)
        {
            rowi_spectrum = map.getSpectrum(rowi);
            // Could handle by sorting, but this shouldn't happen, so want to
            // know if it does
            if (!rowi_spectrum.isSorted()) throw std::runtime_error ("Spectrum not sorted");

            // Iterate over the points in the central spectrum
            for (size_t mzi = 0; mzi < mz_windows; ++mzi)
            {
                // Get the tolerances and value for this point
                double lower_tol_nat = lower_bound_nat[mzi];
                double upper_tol_nat = upper_bound_nat[mzi];
                double lower_tol_iso = lower_bound_iso[mzi];
                double upper_tol_iso = upper_bound_iso[mzi];
                double centre = centre_row_points[mzi].getMZ();
                double sigma = centre * mz_ppm_sigma;
    
                // Select points within tolerance for current spectrum
    
                // Want index of bounds
                Size lower_index = Size(rowi_spectrum.MZBegin(lower_tol_nat) - rowi_spectrum.begin());
                Size upper_index = Size(rowi_spectrum.MZBegin(upper_tol_nat) - rowi_spectrum.begin());
    
                // Check if points found...
                if (lower_index <= upper_index)
                {
                    // Calculate Gaussian value for each found MZ
                    for (Size index = lower_index; index <= upper_index; ++index)
                    {
                        Peak1D peak = rowi_spectrum[index];
                        double mz = peak.getMZ();
                        double intensity = peak.getIntensity();
    
                        // just in case
                        if (mz < lower_tol_nat || mz > upper_tol_nat) continue;
    
                        // calc mz fit
                        mz = (mz - centre) / sigma;
                        mz = -0.5 * mz * mz;
                        double fit = exp(mz) / (sigma * root2pi);
    
                        data_nat[mzi].push_back(intensity);
                        shape_nat[mzi].push_back(fit * rt_shape_nat);
                    }
                }

                // Increment centre to isotope peak
                centre += mz_delta_opt;
                sigma = centre * mz_ppm_sigma;

                // Select points within tolerance for current spectrum
                lower_index = Size(rowi_spectrum.MZBegin(lower_tol_iso) - rowi_spectrum.begin());
                upper_index = Size(rowi_spectrum.MZBegin(upper_tol_iso) - rowi_spectrum.begin());

                // Check if points found...
                if (lower_index <= upper_index)
                {
                    // Calculate Gaussian value for each found MZ
                    for (Size index = lower_index; index <= upper_index; ++index)
                    {
                        Peak1D peak = rowi_spectrum[index];
                        double mz = peak.getMZ();
                        double intensity = peak.getIntensity();

                        // just in case
                        if (mz < lower_tol_iso || mz > upper_tol_iso) continue;

                        // calc mz fit
                        mz = (mz - centre) / sigma;
                        mz = -0.5 * mz * mz;
                        double fit = exp(mz) / (sigma * root2pi);

                        data_iso[mzi].push_back(intensity);
                        shape_iso[mzi].push_back(fit * rt_shape_iso);
                    }
                }
            }
        }
    }

    // Ignore regions with insufficient number of samples
    // If any region is ignored, set all to empty
    for (size_t i = 0; i < mz_windows; ++i) {
        if (data_nat[i].size() < min_sample_opt ||
                       data_iso[i].size() < min_sample_opt) {
            data_nat[i] = {};
            shape_nat[i] = {};
            nAB.push_back(0);
        }
        else
        {
            nAB.push_back(data_nat[i].size() + data_iso[i].size());
        }
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

    // Region means
    double_vect EXa;    // E(Xa) 
    double_vect EXb;    // E(Xb)
    double_vect EYa;    // E(Ya)
    double_vect EYb;    // E(Yb)

    EXa = reduce_2D_vect(data_nat, mean_vector);
    EXb = reduce_2D_vect(data_iso, mean_vector);
    EYa = reduce_2D_vect(shape_nat, mean_vector);
    EYb = reduce_2D_vect(shape_iso, mean_vector);

    // Combined mean is mean of means
    double_vect EXab;    // E(Xab) =  E(E(Xa), E(Xb))
    double_vect EYab;    // E(Yab) =  E(E(Ya), E(Yb))
    EXab = apply_vect_func(EXa, EXb, mean_scalars);
    EYab = apply_vect_func(EYa, EYb, mean_scalars);

    //// FOR TESTING
    // low correlation
    //EXab = EXa;
    //EYab = EYa;
    // high correlation
    //EXab = EXb;
    //EYab = EYb;

    // Centre data in regions relative to combined means
    double_2d CXa;    // (Xa - E(Xab)) --> C(Xa)
    double_2d CXb;    // (Xb - E(Xab))
    double_2d CYa;    // (Ya - E(Yab))
    double_2d CYb;    // (Yb - E(Yab))

    // Centre data in regions relative to combined means
    CXa = apply_vect_func(data_nat, EXab, shift_vector);
    CXb = apply_vect_func(data_iso, EXab, shift_vector);
    CYa = apply_vect_func(shape_nat, EYab, shift_vector);
    CYb = apply_vect_func(shape_iso, EYab, shift_vector);

    // Square vectors
    double_2d CXa2;    // (Xa - E(Xab))^2 --> C(Xa)^2
    double_2d CYa2;    // (Ya - E(Yab))^2

    double_2d CXb2;    // (Xb - E(Xab))^2
    double_2d CYb2;    // (Yb - E(Yab))^2

    CXa2 = apply_vect_func(CXa, square_vector);
    CXb2 = apply_vect_func(CXb, square_vector);
    CYa2 = apply_vect_func(CYa, square_vector);
    CYb2 = apply_vect_func(CYb, square_vector);

    // Products
    double_2d CXaCYa;    // (Xa - E(Xab))(Ya - E(Yab))
    double_2d CXbCYb;    // (Xb - E(Xab))(Yb - E(Yab))

    CXaCYa = apply_vect_func(CXa, CYa, mult_vectors);
    CXbCYb = apply_vect_func(CXb, CYb, mult_vectors);

    // region expected values
    double_vect ECXaCYa;    // E((Xa - E(Xab))(Ya - E(Yab)))
    double_vect ECXa2;    // E((Xa - E(Xab))^2)
    double_vect ECYa2;    // E((Ya - E(Yab))^2)

    double_vect ECXbCYb;    // E((Xb - E(Xab))(Yb - E(Yab)))
    double_vect ECXb2;    // E((Xb - E(Xab))^2)
    double_vect ECYb2;    // E((Yb - E(Yab))^2)

    ECXaCYa = reduce_2D_vect(CXaCYa, mean_vector);
    ECXbCYb = reduce_2D_vect(CXbCYb, mean_vector);
    ECXa2 = reduce_2D_vect(CXa2, mean_vector);
    ECXb2 = reduce_2D_vect(CXb2, mean_vector);
    ECYa2 = reduce_2D_vect(CYa2, mean_vector);
    ECYb2 = reduce_2D_vect(CYb2, mean_vector);

    // Variance, Covariance
    double_vect cov_Xab;
    double_vect var_Xab;
    double_vect var_Yab;

    cov_Xab = apply_vect_func( ECXaCYa, ECXbCYb, mean_scalars);
    var_Xab = apply_vect_func( ECXa2, ECXb2, mean_scalars);
    var_Yab = apply_vect_func( ECYa2, ECYb2, mean_scalars);

    /* Alternate models */
    // low region and high region modelled as flat
    // Region means
    // For high region modelled as all zero, Y_
    // E(Yb) --> E(Y_) = 0 if model region b as all zero

    // Combined mean is mean of means
    double_vect EYa_;    // E(Ya_) =  E(E(Ya), E(Y_)) = 1/2 E(Ya)
    for (const auto& val : EYa) EYa_.push_back(0.5 * val);

    // Centre data in regions relative to combined means
    double_2d CYa_;    // (Ya - E(Ya_))
    // (Yb - E(Ya_)) = -E(Ya_)

    // Centre data in regions relative to combined means
    CYa_ = apply_vect_func(shape_nat, EYa_, shift_vector);

    // Square vectors
    double_2d CYa_2;    // (Ya - E(Ya_))^2
    // (Yb - E(Ya_))^2 = E(Ya_)^2

    CYa_2 = apply_vect_func(CYa_, square_vector);

    // Products
    double_2d CXaCYa_;    // (Xa - E(Xab))(Ya - E(Ya_))
    // (Xb - E(Xab))(Yb - E(Ya_)) --> -E(Ya_)*(Xb - E(Xab))

    CXaCYa_ = apply_vect_func(CXa, CYa_, mult_vectors);

    // region expected values
    double_vect ECXaCYaEa_;    // E((Xa - E(Xab))(Ya - E(Ya_)))
    double_vect ECYaEa_2;    // E((Ya - E(Ya_))^2)

    double_vect ECXbCYbEa_;    // E((Xb - E(Xab))(Yb - E(Ya_))) = -E(Ya_)*E(Xb - E(Xab))
    double_vect ECYbEa_2;    // E((Yb - E(Ya_))^2) = E(Ya_)^2

    ECXaCYaEa_ = reduce_2D_vect(CXaCYa_, mean_vector);
    
    ECXbCYbEa_ = reduce_2D_vect(CXb, mean_vector);    // E(Xb - E(Xab))
    ECXbCYbEa_ = apply_vect_func(ECXbCYbEa_, EYa_, mult_scalars);  // E(Ya_)*E(Xb - E(Xab))
    for (auto& val : ECXbCYbEa_) val = -val;

    ECYaEa_2 = reduce_2D_vect(CYa_2, mean_vector);
    ECYbEa_2 = apply_vect_func(EYa_, EYa_, mult_scalars);  // E(Ya_)^2

    // Variance, Covariance
    double_vect cov_Xa_;
    double_vect var_Ya_;

    cov_Xa_ = apply_vect_func( ECXaCYaEa_, ECXbCYbEa_, mean_scalars);
    var_Ya_ = apply_vect_func( ECYaEa_2, ECYbEa_2, mean_scalars);

    // For low region modelled as all zero, Y_
    // E(Ya) --> E(Y_) = 0 if model region a as all zero

    // Combined mean is mean of means
    double_vect EY_b;    // E(Y_b) =  E(E(Y_), E(Yb)) = 1/2 E(Yb)
    for (const auto& val : EYb) EY_b.push_back(0.5 * val);

    // Centre data in regions relative to combined means
    double_2d CY_b;    // (Yb - E(Y_b))
    // (Ya - E(Y_b)) = -E(Y_b)

    // Centre data in regions relative to combined means
    CY_b = apply_vect_func(shape_iso, EY_b, shift_vector);

    // Square vectors
    double_2d CY_b2;    // (Yb - E(Y_b))^2
    // (Ya - E(Y_b))^2 = E(Y_b)^2

    CY_b2 = apply_vect_func(CY_b, square_vector);

    // Products
    // (Xa - E(Xab))(Ya - E(Y_a)) --> -E(Y_a)*(Xa - E(Xab))
    double_2d CXbCY_b;    // (Xb - E(Xab))(Yb - E(Y_b))

    CXbCY_b = apply_vect_func(CXb, CY_b, mult_vectors);

    // region expected values
    double_vect ECXbCYbE_b;    // E((Xb - E(Xab))(Yb - E(Y_b)))
    double_vect ECYbE_b2;    // E((Yb - E(Y_b))^2)

    double_vect ECXaCYaE_b;    // E((Xa - E(Xab))(Ya - E(Y_b))) = -E(Y_b)*E(Xa - E(Xab))
    double_vect ECYaE_b2;    // E((Ya - E(Y_b))^2) = E(Y_b)^2

    ECXbCYbE_b = reduce_2D_vect(CXbCY_b, mean_vector);
    
    ECXaCYaE_b = reduce_2D_vect(CXa, mean_vector);    // E(Xa - E(Xab))
    ECXaCYaE_b = apply_vect_func(ECXaCYaE_b, EY_b, mult_scalars);  // E(Y_b)*E(Xa - E(Xab))
    for (auto& val : ECXaCYaE_b) val = -val;

    ECYbE_b2 = reduce_2D_vect(CY_b2, mean_vector);
    ECYaE_b2 = apply_vect_func(EY_b, EY_b, mult_scalars);  // E(Y_b)^2

    // Variance, Covariance
    double_vect cov_X_b;
    double_vect var_Y_b;

    cov_X_b = apply_vect_func( ECXaCYaE_b, ECXbCYbE_b, mean_scalars);
    var_Y_b = apply_vect_func( ECYaE_b2, ECYbE_b2, mean_scalars);

    // Correlations
    double_vect correl_XabYab;
    double_vect correl_XabYa_;
    double_vect correl_XabY_b;
    correl_XabYab = correl_vectors(cov_Xab, var_Xab, var_Yab);
    correl_XabYa_ = correl_vectors(cov_Xa_, var_Xab, var_Ya_);
    correl_XabY_b = correl_vectors(cov_X_b, var_Xab, var_Y_b);

    /* correlations between models */
    // Yab, Ya_ covariance...
    // Products
    double_2d CYaCYa_;    // (Ya - E(Yab))(Ya - E(Ya_))
    // (Yb - E(Yab))(Yb - E(Ya_)) --> -E(Ya_)*(Yb - E(Yab))

    CYaCYa_ = apply_vect_func(CYa, CYa_, mult_vectors);

    // region expected values
    double_vect ECYaCYaEa_;    // E((Ya - E(Yab))(Ya - E(Ya_)))
    double_vect ECYbCYbEa_;    // E((Yb - E(Yab))(Yb - E(Ya_))) = -E(Ya_)*E(Yb - E(Yab))

    ECYaCYaEa_ = reduce_2D_vect(CYaCYa_, mean_vector);
    ECYbCYbEa_ = reduce_2D_vect(CYa, mean_vector);    // E(Ya - E(Yab))
    ECYbCYbEa_ = apply_vect_func(ECYbCYbEa_, EYa_, mult_scalars);  // E(Ya_)*E(Ya - E(Yab))
    for (auto& val : ECYbCYbEa_) val = -val;  // -E(Ya_)*E(Ya - E(Yab))

    // Variance, Covariance
    double_vect cov_YabYa_;
    cov_YabYa_ = apply_vect_func( ECYaCYaEa_, ECYbCYbEa_, mean_scalars);

    // Yab, Ya_ correlation
    double_vect correl_YabYa_;
    correl_YabYa_ = correl_vectors(cov_YabYa_, var_Yab, var_Ya_);

    // Yab, Y_b covariance...
    // Products
    double_2d CYbCY_b;    // (Yb - E(Yab))(Yb - E(Y_b))
    // (Ya - E(Yab))(Ya - E(Y_b)) --> -E(Y_b)*(Ya - E(Yab))

    CYbCY_b = apply_vect_func(CYb, CY_b, mult_vectors);

    // region expected values
    double_vect ECYbCYbE_b;    // E((Yb - E(Yab))(Yb - E(Y_b)))
    double_vect ECYaCYaE_b;    // E((Ya - E(Yab))(Ya - E(Y_b))) = -E(Y_b)*E(Ya - E(Yab))

    ECYbCYbE_b = reduce_2D_vect(CYbCY_b, mean_vector);
    ECYaCYaE_b = reduce_2D_vect(CYb, mean_vector);    // E(Yb - E(Yab))
    ECYaCYaE_b = apply_vect_func(ECYaCYaE_b, EY_b, mult_scalars);  // E(Y_b)*E(Yb - E(Yab))
    for (auto& val : ECYaCYaE_b) val = -val;  // -E(Y_b)*E(Yb - E(Yab))

    // Variance, Covariance
    double_vect cov_YabY_b;

    cov_YabY_b = apply_vect_func( ECYaCYaE_b, ECYbCYbE_b, mean_scalars);

    // Yab, Y_b correlation
    double_vect correl_YabY_b;
    correl_YabY_b = correl_vectors(cov_YabY_b, var_Yab, var_Y_b);

    // Compare correlations
    double_vect rm2ABA0;
    double_vect rm2AB0B;
    //double_vect rm2AB1r;

    // Calculate rm values between correlations
    rm2ABA0 = rm_vectors(correl_XabYab, correl_XabYa_);
    rm2AB0B = rm_vectors(correl_XabYab, correl_XabY_b);
    //rm2AB1r = rm_vectors(correlAB, correl1r);

    double_vect fABA0;
    double_vect fAB0B;
    //double_vect fAB1r;

    // Calculate f values between correlation and rm
    fABA0 = f_vectors(correl_XabYa_, rm2ABA0);
    fAB0B = f_vectors(correl_XabY_b, rm2AB0B);
    //fAB1r = f_vectors(correlAB1r, rm2AB1r);

    double_vect hABA0;
    double_vect hAB0B;
    //double_vect hAB1r;

    // Calculate h values between f and rm
    hABA0 = h_vectors(fABA0, rm2ABA0);
    hAB0B = h_vectors(fAB0B, rm2AB0B);
    //hAB1r = h_vectors(fAB1r, rm2AB1r);

    // Subtract 3 and square root
    std::for_each(nAB.begin(), nAB.end(), [](double& d) { d-=3.0;});
    std::transform(nAB.begin(), nAB.end(), nAB.begin(),
                                                 (double(*)(double)) sqrt);

    double_vect zABA0;
    double_vect zAB0B;
    //double_vect zAB1r;

    // Calculate z scores
    zABA0 = z_vectors(correl_XabYab, correl_XabYa_, nAB, correl_YabYa_, hABA0);
    zAB0B = z_vectors(correl_XabYab, correl_XabY_b, nAB, correl_YabY_b, hAB0B);
    //zAB1r = z_vectors(correlAB, correl1r, nAB, correlAB1r, hAB1r);

    double_vect min_score;

    // Find the minimum scores, bounded at zero
    for (size_t idx = 0; idx < zABA0.size(); ++idx) {
        double zA0 = zABA0[idx];
        double z0B = zAB0B[idx];
        //double z1r = zAB1r[idx];
        //double min  = std::min({zA0, zB0, z1r});
        double min  = std::min({zA0, z0B});
        min_score.push_back(std::max({0.0, min}));
    }

    // Package return values
    //double_2d score = {min_score, correlAB, correlA0, correlB0, correl1r};
//    double_2d score = {min_score, {0.0}, {0.0}, {0.0}, {0.0}};

    return min_score;
}
