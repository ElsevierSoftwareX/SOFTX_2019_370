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
   double_2d score;
   for (int n = low_spectra; n < high_spectra; n++)
   {
       score = score_spectra(input_map, n, half_window, opts);

       MSSpectrum<> input_spectrum = input_map.getSpectrum(n);
       MSSpectrum<> output_spectrum = MSSpectrum<Peak1D>(input_spectrum);

       // Copy the computed score into the output intensity
       for (int index = 0; index < input_spectrum.size(); ++index)
       {
           output_spectrum[index].setIntensity(score[0][index]);
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

double_2d
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

    // XXX mz_mu_vect seems like a bad name
    MSSpectrum<> mz_mu_vect = map.getSpectrum(mid_win);

    // Low peak tolerances
    double_vect points_lo_lo;
    double_vect points_lo_hi;
    // High peak tolerances
    double_vect points_hi_lo;
    double_vect points_hi_hi;

    double_2d data_lo;
    double_2d data_hi;
    double_2d shape_lo;
    double_2d shape_hi;
    std::vector<int> len_lo;
    std::vector<int> len_hi;

    // Calculate tolerances for the lo and hi peak for each central MZ
    MSSpectrum<>::Iterator it;
    for (it = mz_mu_vect.begin(); it != mz_mu_vect.end(); ++it)
    {
	double this_mz = it->getMZ();
        points_lo_lo.push_back(this_mz * lo_tol);
        points_lo_hi.push_back(this_mz * hi_tol);
        points_hi_lo.push_back((this_mz + mz_delta_opt) * lo_tol);
        points_hi_hi.push_back((this_mz + mz_delta_opt) * hi_tol);

        double_vect data;
        data_lo.push_back(data);
        data_hi.push_back(data);
        shape_lo.push_back(data);
        shape_hi.push_back(data);
        len_lo.push_back(0);
        len_hi.push_back(0);
    }

    // XXX why float?
    // XXX Can't we compute this once outside the function?
    std::vector<float> rt_shape;

    // Calculate guassian shape in the RT direction
    for (int i = 0; i < (2 * half_window) + 1; ++i)
    {
        float pt = (i - half_window) / rt_sigma;
        pt = -0.5 * pt * pt;
        pt = exp(pt) / (rt_sigma * root2pi);
        rt_shape.push_back(pt);
    }

    // Iterate over the spectra in the window
    for (int rowi = mid_win - half_window; rowi <= mid_win + half_window; ++rowi)
    {
        float rt_lo = rt_shape[rowi - rt_offset];
        float rt_hi = rt_lo;

	MSSpectrum<> rowi_spectrum;
        if (rowi >= 0 && rowi < rt_len)
	{
       	    rowi_spectrum = map.getSpectrum(rowi);
	}

        // Iterate over the points in the central spectrum
        for (size_t mzi = 0; mzi < mz_mu_vect.size(); ++mzi)
	{
            // Get the tolerances and value for this point
            double lo_tol_lo = points_lo_lo[mzi];
            double lo_tol_hi = points_lo_hi[mzi];
            double hi_tol_lo = points_hi_lo[mzi];
            double hi_tol_hi = points_hi_hi[mzi];
            double centre = mz_mu_vect[mzi].getMZ();
            double sigma = centre * mz_ppm_sigma;

            // Check if spectrum within bounds of the file...
            if (rowi >= 0 && rowi < rt_len)
	    {
                // Select points within tolerance for current spectrum
                Size lo_index;
                Size hi_index;
                get_bounds(rowi_spectrum, lo_tol_lo, lo_tol_hi, lo_index, hi_index);

                // Check if points found...
		// XXX should check if the value found at index is near to our target mz
                if (lo_index != -1 && hi_index != -1 && lo_index <= hi_index)
		{
                    // Calculate guassian value for each found MZ
                    for (Size index = lo_index; index <= hi_index; ++index)
		    {
			Peak1D peak = rowi_spectrum[index];
			double mz = peak.getMZ();
			mz = (mz - centre) / sigma;
                        mz = -0.5 * mz * mz;
                        mz = rt_lo * exp(mz) / (sigma * root2pi);
                        shape_lo[mzi].push_back(mz);
		        double intensity = peak.getIntensity();
                        data_lo[mzi].push_back(intensity);
                    }
                    len_lo[mzi] += hi_index - lo_index + 1;

                // ...if not, use dummy data
                }
		else
		{
                    data_lo[mzi].push_back(0);
                    shape_lo[mzi].push_back(rt_lo / (sigma * root2pi));
		    // XXX I think this should be:
		    // len_lo[mzi] += 1;
                }

            // ...if outside use dummy data for this spectrum
            }
	    else
	    {
                data_lo[mzi].push_back(0);
                shape_lo[mzi].push_back(rt_lo / (sigma * root2pi));
		// XXX I think this should be:
		// len_lo[mzi] += 1;
            }

            // Increment centre to hi peak
            centre += mz_delta_opt;
            sigma = centre * mz_ppm_sigma;

            // Check if spectrum within bounds of the file...
            if (rowi >= 0 && rowi < rt_len)
	    {
                // Select points within tolerance for current spectrum
                Size lo_index = -1;
                Size hi_index = -1;
                get_bounds(rowi_spectrum, hi_tol_lo, hi_tol_hi, lo_index, hi_index);

                // Check if points found...
		// XXX should check if the value found at index is near to our target mz
                //if (lo_index <= hi_index)
                if (lo_index != -1 && hi_index != -1 && lo_index <= hi_index)
		{
                    // Calculate guassian value for each found MZ
                    for (Size index = lo_index; index <= hi_index; ++index)
		    {
			Peak1D peak = rowi_spectrum[index];
			double mz = peak.getMZ();
			mz = (mz - centre) / sigma;
                        mz = -0.5 * mz * mz;
                        mz = rt_hi * exp(mz) / (sigma * root2pi);
                        shape_hi[mzi].push_back(mz);
		        double intensity = peak.getIntensity();
                        data_hi[mzi].push_back(intensity);
                    }
                    len_hi[mzi] += hi_index - lo_index + 1;

                // ...if not, use dummy data
                }
		else
		{
                    data_hi[mzi].push_back(0);
                    shape_hi[mzi].push_back(rt_hi / (sigma * root2pi));
		    // XXX I think this should be:
		    // len_hi[mzi] += 1;
                }

            // ...if outside use dummy data for this spectrum
            }
	    else
	    {
                data_hi[mzi].push_back(0);
                shape_hi[mzi].push_back(rt_hi / (sigma * root2pi));
		// XXX I think this should be:
		// len_hi[mzi] += 1;
            }
        }
    }

    //! @todo Combine hi and lo conditions/loops

    // Set point with insufficient samples to zero
    for (size_t leni = 0; leni < len_lo.size(); ++leni) {
        if (len_lo[leni] < min_sample_opt) {
            data_lo[leni]  = {0.0};
            shape_lo[leni] = {0.0};
        }
    }

    // Set points with insufficient samples to zero
    for (size_t leni = 0; leni < len_hi.size(); ++leni) {
        if (len_hi[leni] < min_sample_opt) {
            data_hi[leni]  = {0.0};
            shape_hi[leni] = {0.0};
        } else {
            // Multiply high points by intensity ratio
            for (auto& s : shape_hi[leni]) {
                s *= intensity_ratio_opt;
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
            shape1r_row.push_back(intensity_ratio_opt * length_lo);
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
    double_2d score = {min_score, correlAB, correlA0, correlB0, correl1r};

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
