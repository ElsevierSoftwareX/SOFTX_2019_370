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

#include "pwiz_tools/common/FullReaderList.hpp"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/analysis/spectrum_processing/SpectrumList_MZWindow.hpp"
#include "pwiz/data/msdata/SpectrumInfo.hpp"
#include "pwiz/data/common/cv.hpp"


/*-----------------------------------------------------------------------*/
/******************************* CONSTANTS *******************************/
/*-----------------------------------------------------------------------*/


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
const float default_min_sample      = default_rt_width * default_rt_sigma 
                                        / 2.355;
//! Pi
constexpr double pi() { return std::atan(1) * 4; } 
//! Square root of 2 * Pi
const double root2pi = sqrt(2.0 * pi());


/*-----------------------------------------------------------------------*/
/******************************* TYPEDEFS ********************************/
/*-----------------------------------------------------------------------*/

//! Type definition for a standard vector of doubles
typedef std::vector<double> double_vect;

/*! Type definition for a 2D vector of doubles.
 *
 * Implemented as a standard vector of standard vectors of doubles.
 */
typedef std::vector<double_vect> double_2d;

/*-----------------------------------------------------------------------*/
/******************************** CLASSES ********************************/
/*-----------------------------------------------------------------------*/

/*! @brief Class for holding command line options.
 *
 * This class reads command line arguments and extracts the values. The object
 * can then be passed to all functions that require access to the option.
 */

class Options {

    public:
        const bool getBinaryData = true; //!< Required for pwiz to read data.
        float intensity_ratio; //!< Intensity ratio between lo and hi peaks.
        float rt_width; //!< Retention time FWHM in scans.
        float rt_sigma; //!< Boundary for RT width in SDs.
        float ppm; //!< MZ tolerance in PPM.
        float mz_width; //!< MZ FWHM in PPM.
        float mz_sigma; //!< Boundary for MZ in SDs.
        float mz_delta; //!< MZ difference between peaks.
        float min_sample; //!< Minimum number of points required in each region.
        bool full_out; //!< Output all points (including zero scores).
        std::string mzML_file; //!< Path to input file.
        std::string out_file; //!< Path to output file.

        Options(int argc, char *argv[]);        
};


/*-----------------------------------------------------------------------*/ 
/************************* FUNCTION DECLARATIONS *************************/
/*-----------------------------------------------------------------------*/

//! @brief Print program usage information.
void show_usage(char *cmd);

//! @brief Calculate correlation scores for a window of data.
double_2d score_spectra(pwiz::msdata::MSDataFile &msd, int centre_idx, 
                        int half_window, Options opts);

//! @brief Write correlation scores to an output stream.
void write_scores(double_2d scores, pwiz::msdata::SpectrumPtr raw_data,
                  std::ofstream& out_stream, Options opts); 

//! @brief Centre a vector by subtracting the mean.
double_vect centre_vector(double_vect vect);

//! @brief Square a vector.
double_vect square_vector(double_vect vect);

//! @brief Sum a vector.
double sum_vector(double_vect vect);

//! @brief Elementwise multiplication of two vectors.
double_vect mult_vectors(double_vect vect1, double_vect vect2);

//! @brief Elementwise division of two vectors.
double_vect div_vectors(double_vect vect1, double_vect vect2);

//! @brief Calculate correlation between vectors.
double_vect correl_vectors(double_vect vect1, double_vect vect2, 
                           double_vect vect3);

//! @brief Calculate rm values from two vectors.
double_vect rm_vectors(double_vect vect1, double_vect vect2);

//! @brief Calculate f values from two vectors.
double_vect f_vectors(double_vect correl_vect, double_vect rm_vect);

//! @brief Calculate h values from two vectors.
double_vect h_vectors(double_vect f_vect, double_vect rm_vect);

//! @brief Calculate z values.
double_vect z_vectors(double_vect cor1, double_vect cor2, double_vect sqrtn,
                      double_vect cross_cor, double_vect h_vect);

//! @brief Apply function to each element of a vector. 
template <typename T, typename F>
std::vector<T> apply_vect_func(std::vector<T> vect, F func);

//! @brief Apply function that combines elements of two vectors.
template <typename T, typename F>
std::vector<T> apply_vect_func(std::vector<T> vect1, std::vector<T> vect2,
                                                                    F func);
//! @brief Apply function that reduces a 2D vector to a 1D vector.
template <typename T, typename F>
std::vector<T> reduce_2D_vect (std::vector<std::vector<T>> vect2D, F func);

 
/*-----------------------------------------------------------------------*/
/********************************* MAIN **********************************/
/*-----------------------------------------------------------------------*/


int main(int argc, char *argv[])
{
    // Read user options
    Options opts(argc, argv);

    // Calculate number of spectra for each window
    int half_window = ceil(opts.rt_sigma * opts.rt_width / 2.355);

    // Access the input file using pwiz
    pwiz::msdata::FullReaderList readers;
    pwiz::msdata::MSDataFile msd(opts.mzML_file, &readers);
    pwiz::msdata::SpectrumList& spectrumList = *msd.run.spectrumListPtr;

    int rt_len = spectrumList.size();
    double_2d score;
    pwiz::msdata::SpectrumPtr centre_vect;
    
    // Setup output stream
    std::ofstream outfile;
    outfile.open(opts.out_file);
    outfile.precision(12);

    // Loop over the spectra
    for (int centre_rt = 0; centre_rt < rt_len; ++centre_rt) {

        // Select the centre spectrum
        centre_vect = spectrumList.spectrum(centre_rt, opts.getBinaryData);
        
        // Score the window
        score = score_spectra(msd, centre_rt, half_window, opts);
        
        // Output the results
        write_scores(score, centre_vect, outfile, opts);
    }

    // Close ouput stream
    outfile.close();

    std::cout << "Done!" << std::endl;
    return 0;
}

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


/*-----------------------------------------------------------------------*/
/************************* FUNCTION DEFINITIONS **************************/
/*-----------------------------------------------------------------------*/

/*! Display program usage information to the user. Called when no arguments
 * are supplied or the -h option is given.
 *
 * @param cmd The original command given by the user.
 */
void show_usage(char *cmd)
{
    using namespace std;

    cout << "Usage:     " << cmd << " [-options] [arguments]"       << endl;
    cout                                                            << endl;
    cout << "options:   " << "-h  show this help information"       << endl;
    cout << "           " << "-i  ratio of doublet intensities (isotope \n";
    cout << "           " << "    / parent)"                        << endl;
    cout << "           " << "-r  full width at half maximum for \n"       ;
    cout << "           " << "    retention time in number of scans"<< endl;
    cout << "           " << "-R  retention time width boundary in \n"     ;
    cout << "           " << "    standard deviations"              << endl;
    cout << "           " << "-p  m/z tolerance in parts per million"      ;
    cout                                                            << endl;
    cout << "           " << "-m  m/z full width at half maximum in \n"    ;
    cout << "           " << "    parts per million"                << endl;
    cout << "           " << "-M  m/z window boundary in standard \n"      ;
    cout << "           " << "    deviations"                       << endl;
    cout << "           " << "-D  m/z difference for doublets"      << endl;
    cout << "           " << "-s  minimum number of data points \n"        ;
    cout << "           " << "    required in each sample region"   << endl;
    cout << "           " << "-o  turn on full output, including zero \n"  ;
    cout << "           " << "    score points"                     << endl;
    cout                                                            << endl;
    cout << "arguments: " << "mzML_file     path to mzML file"      << endl;
    cout << "           " << "out_file      path to output file"    << endl;
    cout                                                            << endl;
    cout << "example:   " << cmd << " example.mzML output.txt"      << endl;
    cout                                                            << endl;
}

/*! Scores are output in CSV format with the following fields: retention time,
 * mz, intensity, minimum score, correlAB, correlA0, correlB0, correl1r.
 *
 * @param scores 2D vector of scores returned by score_spectra.
 * @param raw_data Spectrum pointer to the raw central vector.
 * @param out_stream Stream to write output too.
 * @param opts User defined Options object.
 */
void write_scores(double_2d scores, pwiz::msdata::SpectrumPtr raw_data,
                  std::ofstream& out_stream, Options opts)
{
    // Get central spectrum retention time
    pwiz::msdata::SpectrumInfo spectrum_info;
    spectrum_info.update(*raw_data, opts.getBinaryData);
    double rt = spectrum_info.retentionTime;
    
    // Get raw MZ/Intensity pairs
    std::vector<pwiz::msdata::MZIntensityPair> raw_pairs;
    raw_data->getMZIntensityPairs(raw_pairs);
    
    // Write output
    for (size_t idx = 0; idx < raw_pairs.size(); ++idx) {
        double mz  = raw_pairs[idx].mz;
        double amp = raw_pairs[idx].intensity;
        double ms  = scores[0][idx]; 
        double AB  = scores[1][idx];
        double A0  = scores[2][idx];
        double B0  = scores[3][idx];
        double r1  = scores[4][idx];

        if (opts.full_out == true) {
            out_stream << rt << ", " << mz << ", " << amp << ", " 
                       << ms << ", " << AB << ", " << A0 << ", " 
                       << B0 << ", " << r1 << std::endl; 
        } else {
            if (ms > 0.0) {
                out_stream << rt << ", " << mz << ", " << amp << ", " 
                           << ms << ", " << AB << ", " << A0 << ", " 
                           << B0 << ", " << r1 << std::endl;
            }
        }
    }
}

/*! Calculate the mean of all values in the vector and subtract from each
 * individual value.
 *
 * @param vect The vector to centre.
 *
 * @return The centred vector.
 */
double_vect centre_vector(double_vect vect)
{
    double sum  = std::accumulate(vect.begin(), vect.end(), 0.0);
    double mean = sum / vect.size();
    double_vect centered;

    for (auto v : vect) {
        centered.push_back(v - mean);
    }

    return centered;
}

/*! Multiply each value in a vector by itself.
 *
 * @param vect The vector to square.
 *
 * @return The squared vector.
 */
double_vect square_vector(double_vect vect)
{
    double_vect squared;

    for (auto v : vect) {
        squared.push_back(v * v);
    }

    return squared;
}

/*! Add up all the values in a vector.
 *
 * @param vect The vector to sum. 
 *
 * @return The summed value. 
 */
double sum_vector(double_vect vect)
{
    double sum = std::accumulate(vect.begin(), vect.end(), 0.0);

    return sum;
}

/*! Multiply pairs of values from two vectors of equal length.
 *
 * @param vect1 First vector to be multiplied.
 * @param vect2 Second vector to be multiplied.
 *
 * @return Vector of pairwise multiples.
 */
double_vect mult_vectors(double_vect vect1, double_vect vect2)
{
    double_vect mult;

    // Throw exception if vectors of different lengths
    if (vect1.size() != vect2.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    // Pairwise multiplication
    for (size_t idx = 0; idx < vect1.size(); ++idx) {
        mult.push_back(vect1[idx] * vect2[idx]);
    }

    return mult;
}
   
/*! Divide values from one vector by values from a second vector of equal 
 * length.
 *
 * @param vect1 The vector to be divided (numerator).
 * @param vect2 The vector to divide by (denominator).
 *
 * @return Vector of vect1[i] / vect2[i] values.
 */
double_vect div_vectors(double_vect vect1, double_vect vect2)
{
    double_vect divided;

    // Throw exception if vectors of different lengths
    if (vect1.size() != vect2.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    // Pairwise division
    for (size_t idx = 0; idx < vect1.size(); ++idx) {
        divided.push_back(vect1[idx] / vect2[idx]);
    }

    return divided;

}

/*!@todo Explain input vectors and return correl vector
 */
double_vect correl_vectors(double_vect vect1, double_vect vect2, 
                           double_vect vect3)
{
    double_vect correlated;
    double_vect mult;

    mult = mult_vectors(vect2, vect3);
    std::transform(mult.begin(), mult.end(), mult.begin(), 
                                            (double(*)(double)) std::sqrt);
    correlated = div_vectors(vect1, mult);

    for (auto& c : correlated) {
        if(isnan(c)) {
            c = 0;
        }
        if(c < 0) {
            c = 0;
        }
    }

    return correlated;
}

/*! Calculate values equal to 0.5 * (vect1[i]^2 + vect2[i]^2)
 *
 * @param vect1 The first values in the equation.
 * @param vect2 The second values in the equation.
 *
 * @return Vector of calculated RM values. 
 */
double_vect rm_vectors(double_vect vect1, double_vect vect2)
{
    double_vect rm;

    // Throw exception if vectors of different lengths
    if (vect1.size() != vect2.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    // Calculate RM values.
    for (size_t idx = 0; idx < vect1.size(); ++idx) {
        double tmp = (vect1[idx] * vect1[idx]) + (vect2[idx] * vect2[idx]);
        rm.push_back(0.5 * tmp);
    }

    return rm;

}

/*!F values are bounded at an upper value of 1.0.
 *
 * @param correl_vect Vector returned by correl_vectors.
 * @param rm_vect Vector return by rm_vectors.
 *
 * @return Vector of calculated f values.
 *
 * @todo Explain what f value is.
 */
double_vect f_vectors(double_vect correl_vect, double_vect rm_vect)
{
    double_vect f_vect;

    // Throw exception if vectors of different lengths
    if (correl_vect.size() != rm_vect.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    // Calculate f values
    for (size_t idx = 0; idx < correl_vect.size(); ++idx) {
        double correl = correl_vect[idx];
        double rm     = rm_vect[idx];

        f_vect.push_back((1.0 - correl) / (2.0 * (1.0 - rm)));
    }

    // Set upper bound to 1.0
    for (auto& f : f_vect) {
        if(f > 1.0) {
            f = 1.0;
        }
    }

    return f_vect;
}

/*! @param f_vect Vector of values returned by f_vectors. 
 * @param rm_vect Vector of values returned by rm_vectors.
 *
 * @return Vector of calculated h values.
 * 
 * @todo Explain what h value is.
 */
double_vect h_vectors(double_vect f_vect, double_vect rm_vect)
{
    double_vect h_vect;

    // Throw exception if vectors of different lengths
    if (f_vect.size() != rm_vect.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    // Calculate h value
    for (size_t idx = 0; idx < f_vect.size(); ++idx) {
        double f  = f_vect[idx];
        double rm = rm_vect[idx];

        h_vect.push_back((1.0 - f * rm) / (1.0 - rm));
    }

    return h_vect;
}

/*! @todo Explain z score calculation and params
 */
double_vect z_vectors(double_vect cor1, double_vect cor2, double_vect sqrtn,
                      double_vect cross_cor, double_vect h_vect)
{
    double_vect z_vect;

    // Throw exception if vectors of different lengths
    if (cor1.size() != cor2.size()      ||
        cor1.size() != sqrtn.size()     ||
        cor1.size() != cross_cor.size() ||
        cor1.size() != h_vect.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }
    
    // Calculate z score
    for (size_t idx = 0; idx < cor1.size(); ++idx) {
        
        double z1  = std::atanh(cor1[idx]);
        double z2  = std::atanh(cor2[idx]);
        
        double num   = (z1 - z2) * sqrtn[idx];
        double denom = 2.0 * (1.0 - cross_cor[idx]) * h_vect[idx];

        z_vect.push_back(num / std::sqrt(denom));
    }

    return z_vect;
}

/*! Template function that takes a vector of any type and applies a function
 * to each element to give another vector of the same type.
 *
 * @param vect Vector to apply the function too.
 * @param func Function to apply to each element of the vector. Must return a
 * value of same type as the vector.
 *
 * @returns Vector containing results from applying _func_ to _vect_.
 */
template <typename T, typename F>
std::vector<T> apply_vect_func(std::vector<T> vect, F func)
{
    std::vector<T> applied;
    
    for (auto v : vect) {
        applied.push_back(func(v));
    }

    return applied;
}

/*! Template function that takes two vectors of any type and applies a function
 * to each pair of elements to give another vector of the same type.
 *
 * @param vect1 First vector to apply the function too.
 * @param vect2 Second vector to apply the function too
 * @param func Function to apply to each pair of elements. Must return a single
 * value of the same type as the vectors
 *
 * @returns Vector containing results from applying _func_ to _vect1_ and 
 * _vect2_.
 */
template <typename T, typename F>
std::vector<T> apply_vect_func(std::vector<T> vect1, std::vector<T> vect2, 
                                                                     F func)
{
    std::vector<T> applied;
    
    // Throw exception if vectors of different lengths
    if (vect1.size() != vect2.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    for (size_t idx = 0; idx < vect1.size(); ++idx) {
        applied.push_back(func(vect1[idx], vect2[idx]));
    }

    return applied;
}

/*! Template function that takes a 2D vector of any type and applies a function
 * to each element vector to give a 1D vector of the same type.
 *
 * @param vect2D 2D vector to apply the function too.
 * @param func Function to apply to each element of the 2D vector. Must take a 
 * vector and return a single value of same type as the vector.
 *
 * @returns Vector containing results from applying _func_ to _vect2D_.
 */
template <typename T, typename F>
std::vector<T> reduce_2D_vect (std::vector<std::vector<T>> vect2D, F func)
{
    std::vector<T> reduced;

    for (auto vect : vect2D) {
        reduced.push_back(func(vect));
    }

    return reduced;
}

/*-----------------------------------------------------------------------*/
/***************************** CLASS METHODS *****************************/
/*-----------------------------------------------------------------------*/

/*! @brief Options object constructor.
 *
 * Construct new Options object by reading in arguments from the command line.
 * The number of required arguments in checked.
 * 
 * @param argc The argument count passed to Main
 * @param argv The argument value array passed to Main
 *
 * @todo Validate user input option values
 */

Options::Options(int argc, char *argv[])
{
    char opt;
    int opt_idx;

    intensity_ratio = default_intensity_ratio;
    rt_width        = default_rt_width;
    rt_sigma        = default_rt_sigma;
    ppm             = default_ppm;
    mz_width        = default_fwhm;
    mz_sigma        = default_mz_sigma;
    mz_delta        = default_mz_delta;
    min_sample      = default_min_sample;
    full_out        = false;
    mzML_file       = "";
    out_file        = "";

    // Show usage and exit if no options are given
    if (argc == 1) {
        show_usage(argv[0]);
        exit(1);
    }

    // Check arguments and assign to attributes
    while ((opt = getopt(argc, argv, "hd:i:r:R:p:m:M:D:s:o")) != -1){
        
        switch (opt) {
            case 'h':
                show_usage(argv[0]);
                exit(1);
                break;
            case 'i':
                intensity_ratio = std::stof(std::string(optarg));
                break;
            case 'r':
                rt_width = std::stof(std::string(optarg));
                break;
            case 'R':
                rt_sigma = std::stof(std::string(optarg));
                break;
            case 'p':
                ppm = std::stof(std::string(optarg));
                break;
            case 'm':
                mz_width = std::stof(std::string(optarg));
                break;
            case 'M':
                mz_sigma = std::stof(std::string(optarg));
                break;
            case 'D':
                mz_delta = std::stof(std::string(optarg));
                break;
            case 's':
                min_sample = std::stof(std::string(optarg));
                break;
            case 'o':
                full_out = true;
                break;
        }
    }

    // Read remaining text arguments
    for (opt_idx = optind; opt_idx < argc; opt_idx++) {

        if (mzML_file == "") { 
            mzML_file = argv[opt_idx];
        } else if (out_file == "") {
            out_file = argv[opt_idx];
        } else {
            std::cout << "Too many arguments supplied. See usage.";
            std::cout << std::endl;
            exit(1);
        }
    }

    // Check that all attributes have been set
    if (out_file == "") {
        std::cout << "Insufficient arguments supplies. See usage.";
        std::cout << std::endl;
        exit(1);
    }
}
