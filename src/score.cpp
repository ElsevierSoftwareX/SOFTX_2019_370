/*
#include <iostream>
#include <unistd.h>
#include <math.h>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <functional>
#include <fstream>
#include <numeric>
#include <iostream>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include "vector.h"

using namespace OpenMS;
using namespace std;
*/


//! @brief Calculate correlation scores for a window of data.

/*
double_2d score_spectra(pwiz::msdata::MSDataFile &msd, int centre_idx,
                        int half_window, Options opts);

//! @brief Write correlation scores to an output stream.
void write_scores(double_2d scores, pwiz::msdata::SpectrumPtr raw_data,
                  std::ofstream& out_stream, Options opts);
*/


/*! Calculate correlation scores for each MZ point in a central spectrum of
 * a data window.
 *
 * @param msd Reference to the input MS data file.
 * @param centre_idx Index of the spectrum to score.
 * @param half_window The number of spectra each side of the central spectrum
 * to include.
 * @param opts Options object.
 *
 * @return Vector of five vectors (min score, correlAB, correlA0, correlB0,
 * correl1r) giving score for each MZ in the central spectrum.
 */

/*
double_2d
score_spectra(pwiz::msdata::MSDataFile &msd, int centre_idx,
              int half_window, Options opts)
{
    // Create pointer to the spectrum list
    pwiz::msdata::SpectrumList& spectrumList = *msd.run.spectrumListPtr;

    // Calculate constant values
    float  rt_sigma     = opts.rt_width / 2.355;
    double mz_ppm_sigma = opts.mz_width / 2.355e6;
    int    rt_len       = spectrumList.size();
    int    mid_win      = centre_idx;
    double lo_tol       = 1.0 - opts.mz_sigma * mz_ppm_sigma;
    double hi_tol       = 1.0 + opts.mz_sigma * mz_ppm_sigma;
    int    rt_offset    = mid_win - half_window;

    // Extract the central spectrum
    std::vector<pwiz::msdata::MZIntensityPair> mz_mu_pairs;
    pwiz::msdata::SpectrumPtr mz_mu_vect;
    mz_mu_vect = spectrumList.spectrum(mid_win, opts.getBinaryData);
    mz_mu_vect->getMZIntensityPairs(mz_mu_pairs);

    double_vect points_lo_lo;
    double_vect points_lo_hi;
    double_vect points_hi_lo;
    double_vect points_hi_hi;

    double_2d data_lo;
    double_2d data_hi;
    double_2d shape_lo;
    double_2d shape_hi;
    std::vector<int> len_lo;
    std::vector<int> len_hi;

    // Calculate tolerances for the lo and hi peak for each central MZ
    for (auto pair : mz_mu_pairs) {
        points_lo_lo.push_back(pair.mz * lo_tol);
        points_lo_hi.push_back(pair.mz * hi_tol);
        points_hi_lo.push_back((pair.mz + opts.mz_delta) * lo_tol);
        points_hi_hi.push_back((pair.mz + opts.mz_delta) * hi_tol);

        double_vect data;
        data_lo.push_back(data);
        data_hi.push_back(data);
        shape_lo.push_back(data);
        shape_hi.push_back(data);
        len_lo.push_back(0);
        len_hi.push_back(0);
    }

    std::vector<float> rt_shape;

    // Calculate guassian shape in the RT direction
    for (int i = 0; i < (2 * half_window) + 1; ++i) {

        float pt = (i - half_window) / rt_sigma;
        pt = -0.5 * pt * pt;
        pt = exp(pt) / (rt_sigma * root2pi);

        rt_shape.push_back(pt);
    }

    // Iterate over the spectra in the window
    for (int rowi = mid_win - half_window;
         rowi <= mid_win + half_window; ++rowi) {

        float rt_lo = rt_shape[rowi - rt_offset];
        float rt_hi = rt_lo;

        // Iterate over the points in the central spectrum
        for (size_t mzi = 0; mzi < mz_mu_pairs.size(); ++mzi) {

            // Get the tolerances and value for this point
            double lo_tol_lo = points_lo_lo[mzi];
            double lo_tol_hi = points_lo_hi[mzi];
            double hi_tol_lo = points_hi_lo[mzi];
            double hi_tol_hi = points_hi_hi[mzi];
            double centre    = mz_mu_pairs[mzi].mz;
            double sigma     = centre * mz_ppm_sigma;

            // Use pwiz to get point with tolerance from spectra
            pwiz::analysis::SpectrumList_MZWindow lo_window(
                                                msd.run.spectrumListPtr,
                                                lo_tol_lo, lo_tol_hi);
            pwiz::analysis::SpectrumList_MZWindow hi_window(
                                                msd.run.spectrumListPtr,
                                                hi_tol_lo, hi_tol_hi);
            pwiz::msdata::SpectrumPtr lo_spectrum;
            pwiz::msdata::SpectrumPtr hi_spectrum;
            std::vector<pwiz::msdata::MZIntensityPair> lo_pairs;
            std::vector<pwiz::msdata::MZIntensityPair> hi_pairs;

            // Check if spectrum within bounds of the file...
            if (rowi >= 0 && rowi < rt_len) {

                // Select points within tolerance for current spectrum
                lo_spectrum = lo_window.spectrum(rowi, opts.getBinaryData);
                lo_spectrum->getMZIntensityPairs(lo_pairs);

                // Check if points found...
                if (lo_pairs.size() > 0) {

                    // Calculate guassian value for each found MZ
                    for (auto pair : lo_pairs) {
                        float mz = (pair.mz - centre) / sigma;
                        mz = -0.5 * mz * mz;
                        mz = rt_lo * exp(mz) / (sigma * root2pi);
                        shape_lo[mzi].push_back(mz);
                        data_lo[mzi].push_back(pair.intensity);
                    }

                    len_lo[mzi] += lo_pairs.size();

                // ...if not, use dummy data
                } else {
                    data_lo[mzi].push_back(0);
                    shape_lo[mzi].push_back(rt_lo / (sigma * root2pi));
                }

            // ...if outside use dummy data for this spectrum
            } else {
                data_lo[mzi].push_back(0);
                shape_lo[mzi].push_back(rt_lo / (sigma * root2pi));
            }

            // Increment centre to hi peak
            centre += opts.mz_delta;
            sigma = centre * mz_ppm_sigma;

            // Check if spectrum within bounds of the file...
            if (rowi >= 0 && rowi < rt_len) {

                // Select points within tolerance for current spectrum
                hi_spectrum = hi_window.spectrum(rowi, opts.getBinaryData);
                hi_spectrum->getMZIntensityPairs(hi_pairs);

                // Check if points found...
                if (hi_pairs.size() > 0) {

                    // Calculate guassian value for each found MZ
                    for (auto pair : hi_pairs) {
                        float mz = (pair.mz - centre) / sigma;
                        mz = -0.5 * mz * mz;
                        mz = rt_hi * exp(mz) / (sigma * root2pi);
                        shape_hi[mzi].push_back(mz);
                        data_hi[mzi].push_back(pair.intensity);
                    }
                    len_hi[mzi] += hi_pairs.size();

                // ...if not, use dummy data
                } else {
                    data_hi[mzi].push_back(0);
                    shape_hi[mzi].push_back(rt_hi / (sigma * root2pi));
                }

            // ...if outside use dummy data for this spectrum
            } else {
                data_lo[mzi].push_back(0);
                shape_lo[mzi].push_back(rt_hi / (sigma * root2pi));
            }
        }
    }

    //! @todo Combine hi and lo conditions/loops

    // Set point with insufficient samples to zero
    for (size_t leni = 0; leni < len_lo.size(); ++leni) {
        if (len_lo[leni] < opts.min_sample) {
            data_lo[leni]  = {0.0};
            shape_lo[leni] = {0.0};
        }
    }

    // Set points with insufficient samples to zero
    for (size_t leni = 0; leni < len_hi.size(); ++leni) {
        if (len_hi[leni] < opts.min_sample) {
            data_hi[leni]  = {0.0};
            shape_hi[leni] = {0.0};
        } else {
            // Multiply high points by intensity ratio
            for (auto& s : shape_hi[leni]) {
                s *= opts.intensity_ratio;
            }
        }
    }

    double_2d dataAB;
    double_vect nAB;

    // Setup combined data vector
    for (size_t i = 0; i < data_lo.size(); ++i) {

        double_vect dataAB_row;
        size_t length_lo = data_lo[i].size();
        size_t length_hi = data_hi[i].size();

        for (auto lo_value : data_lo[i]){
            dataAB_row.push_back(lo_value * length_hi);
        }
        for (auto hi_value : data_hi[i]){
            dataAB_row.push_back(hi_value * length_lo);
        }
        dataAB.push_back(dataAB_row);

        nAB.push_back(length_lo + length_hi);
    }

    double_2d shapeAB;
    double_2d shapeA0;
    double_2d shapeB0;
    double_2d shape1r;

    // Setup alternative shape vectors
    for (size_t i = 0; i < shape_lo.size(); ++i) {

        double_vect shapeAB_row;
        double_vect shapeA0_row;
        double_vect shapeB0_row;
        double_vect shape1r_row;
        size_t length_lo = shape_lo[i].size();
        size_t length_hi = shape_hi[i].size();

        for (auto lo_value : shape_lo[i]){
            shapeAB_row.push_back(lo_value * length_hi);
            shapeA0_row.push_back(lo_value * length_hi);
            shapeB0_row.push_back(0.0);
            shape1r_row.push_back(length_hi);
        }
        for (auto hi_value : shape_hi[i]){
            shapeAB_row.push_back(hi_value * length_lo);
            shapeA0_row.push_back(0.0);
            shapeB0_row.push_back(hi_value * length_lo);
            shape1r_row.push_back(opts.intensity_ratio * length_lo);
        }

        shapeAB.push_back(shapeAB_row);
        shapeA0.push_back(shapeA0_row);
        shapeB0.push_back(shapeB0_row);
        shape1r.push_back(shape1r_row);
    }

    // Centre vectors
    dataAB  = apply_vect_func(dataAB,  centre_vector);
    shapeAB = apply_vect_func(shapeAB, centre_vector);
    shapeA0 = apply_vect_func(shapeA0, centre_vector);
    shapeB0 = apply_vect_func(shapeB0, centre_vector);
    shape1r = apply_vect_func(shape1r, centre_vector);

    double_2d data2AB;
    double_2d shape2AB;
    double_2d shape2A0;
    double_2d shape2B0;
    double_2d shape21r;

    // Square vectors
    data2AB  = apply_vect_func(dataAB,  square_vector);
    shape2AB = apply_vect_func(shapeAB, square_vector);
    shape2A0 = apply_vect_func(shapeA0, square_vector);
    shape2B0 = apply_vect_func(shapeB0, square_vector);
    shape21r = apply_vect_func(shape1r, square_vector);

    double_vect SSY;
    double_vect SSXAB;
    double_vect SSXA0;
    double_vect SSXB0;
    double_vect SSX1r;

    // Sum squared vectors
    SSY   = reduce_2D_vect(data2AB,  sum_vector);
    SSXAB = reduce_2D_vect(shape2AB, sum_vector);
    SSXA0 = reduce_2D_vect(shape2A0, sum_vector);
    SSXB0 = reduce_2D_vect(shape2B0, sum_vector);
    SSX1r = reduce_2D_vect(shape21r, sum_vector);

    double_2d datashape;
    double_vect SXYAB;
    double_vect SXYA0;
    double_vect SXYB0;
    double_vect SXY1r;
    double_vect SXYABA0;
    double_vect SXYABB0;
    double_vect SXYAB1r;

    // Sum cross products of alternative shapes
    datashape = apply_vect_func(dataAB, shapeAB, mult_vectors);
    SXYAB     = reduce_2D_vect(datashape, sum_vector);
    datashape = apply_vect_func(dataAB, shapeA0, mult_vectors);
    SXYA0     = reduce_2D_vect(datashape, sum_vector);
    datashape = apply_vect_func(dataAB, shapeB0, mult_vectors);
    SXYB0     = reduce_2D_vect(datashape, sum_vector);
    datashape = apply_vect_func(dataAB, shape1r, mult_vectors);
    SXY1r     = reduce_2D_vect(datashape, sum_vector);
    datashape = apply_vect_func(shapeAB, shapeA0, mult_vectors);
    SXYABA0   = reduce_2D_vect(datashape, sum_vector);
    datashape = apply_vect_func(shapeAB, shapeB0, mult_vectors);
    SXYABB0   = reduce_2D_vect(datashape, sum_vector);
    datashape = apply_vect_func(shapeAB, shape1r, mult_vectors);
    SXYAB1r   = reduce_2D_vect(datashape, sum_vector);

    double_vect correlAB;
    double_vect correlA0;
    double_vect correlB0;
    double_vect correl1r;
    double_vect correlABA0;
    double_vect correlABB0;
    double_vect correlAB1r;

    // Calcualate correlation between shapes
    correlAB   = correl_vectors(SXYAB,   SSXAB, SSY);
    correlA0   = correl_vectors(SXYA0,   SSXA0, SSY);
    correlB0   = correl_vectors(SXYB0,   SSXB0, SSY);
    correl1r   = correl_vectors(SXY1r,   SSX1r, SSY);
    correlABA0 = correl_vectors(SXYABA0, SSXAB, SSXA0);
    correlABB0 = correl_vectors(SXYABB0, SSXAB, SSXB0);
    correlAB1r = correl_vectors(SXYAB1r, SSXAB, SSX1r);

    double_vect rm2ABA0;
    double_vect rm2ABB0;
    double_vect rm2AB1r;

    // Calculate rm values between correlations
    rm2ABA0 = rm_vectors(correlAB, correlA0);
    rm2ABB0 = rm_vectors(correlAB, correlB0);
    rm2AB1r = rm_vectors(correlAB, correl1r);

    double_vect fABA0;
    double_vect fABB0;
    double_vect fAB1r;

    // Calculate f values between correlation and rm
    fABA0 = f_vectors(correlABA0, rm2ABA0);
    fABB0 = f_vectors(correlABB0, rm2ABB0);
    fAB1r = f_vectors(correlAB1r, rm2AB1r);

    double_vect hABA0;
    double_vect hABB0;
    double_vect hAB1r;

    // Calculate h values between f and rm
    hABA0 = h_vectors(fABA0, rm2ABA0);
    hABB0 = h_vectors(fABB0, rm2ABB0);
    hAB1r = h_vectors(fAB1r, rm2AB1r);

    // Subtract 3 and square root
    std::for_each(nAB.begin(), nAB.end(), [](double& d) { d-=3.0;});
    std::transform(nAB.begin(), nAB.end(), nAB.begin(),
                                                 (double(*)(double)) sqrt);

    double_vect zABA0;
    double_vect zABB0;
    double_vect zAB1r;

    // Calculate z scores
    zABA0 = z_vectors(correlAB, correlA0, nAB, correlABA0, hABA0);
    zABB0 = z_vectors(correlAB, correlB0, nAB, correlABB0, hABB0);
    zAB1r = z_vectors(correlAB, correl1r, nAB, correlAB1r, hAB1r);

    double_vect min_score;

    // Find the minimum scores, bounded at zero
    for (size_t idx = 0; idx < zABA0.size(); ++idx) {
        double zA0 = zABA0[idx];
        double zB0 = zABB0[idx];
        double z1r = zAB1r[idx];
        double min  = std::min({zA0, zB0, z1r});
        min_score.push_back(std::max({0.0, min}));
    }

    // Package return values
    double_2d score = {min_score, correlAB, correlA0,
                                              correlB0, correl1r};

    return score;
}
*/
