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


/*-----------------------------------------------------------------------*/
/******************************* CONSTANTS *******************************/
/*-----------------------------------------------------------------------*/


// default difference in mass of isotopes
const float default_mz_delta        = 6.0201;
// default m/z tolerance in parts per million
const float default_ppm             = 4.0;
// Full Width Half Maximum in PPM
const float default_fwhm            = 150.0;
const float default_mz_sigma        = 1.5;
// default ratio of isotopes
const float default_intensity_ratio = 1.0;
// default retention time FWHM in scans 
const float default_rt_width        = 17.0;
const float default_rt_sigma        = 1.5;
// minimum number of samples in score regions
const float default_min_sample      = default_rt_width * default_rt_sigma 
                                        / 2.355;
// pi
constexpr double pi() { return std::atan(1) * 4; } 
// sqrt 2pi
const double root2pi = sqrt(2.0 * pi());

/*-----------------------------------------------------------------------*/ 
/************************* FUNCTION DECLARATIONS *************************/
/*-----------------------------------------------------------------------*/

void show_usage(char *cmd);

std::vector<double> centre_vector(std::vector<double> vect);

std::vector<double> square_vector(std::vector<double> vect);

double sum_vector(std::vector<double> vect);

std::vector<double> mult_vectors(std::vector<double> vect1, 
                                                std::vector<double> vect2);

std::vector<double> div_vectors(std::vector<double> vect1, 
                                                std::vector<double> vect2);
  
std::vector<double> correl_vectors(std::vector<double> vect1,
                     std::vector<double> vect2, std::vector<double> vect3);

std::vector<double> rm_vectors(std::vector<double> vect1, 
                                                std::vector<double> vect2);

std::vector<double> f_vectors(std::vector<double> correl_vect,
                                              std::vector<double> rm_vect);

std::vector<double> h_vectors(std::vector<double> f_vect,
                                              std::vector<double> rm_vect);

std::vector<double> z_vectors(std::vector<double> cor1,
                              std::vector<double> cor2,
                              std::vector<double> sqrtn,
                              std::vector<double> cross_cor,
                              std::vector<double> h_vect);

template <typename T, typename F>
std::vector<T> apply_vect_func(std::vector<T> vect, F func);

template <typename T, typename F>
std::vector<T> apply_vect_func(std::vector<T> vect1, std::vector<T> vect2,
                                                                    F func);

//template <typename T, typename F>
//std::vector<std::vector<T>> apply_vect_func(
//                                  std::vector<std::vector<T>> vect, F func);

template <typename T, typename F>
std::vector<T> reduce_2D_vect (std::vector<std::vector<T>> vect2D, F func);

  
/*-----------------------------------------------------------------------*/
/******************************** CLASSES ********************************/
/*-----------------------------------------------------------------------*/


class Options {

    public:
        float intensity_ratio;
        float rt_width;
        float rt_sigma;
        float ppm;
        float mz_width;
        float mz_sigma;
        float mz_delta;
        float min_sample;
        std::string mzML_file;

        Options(int argc, char *argv[]);        
};


/*-----------------------------------------------------------------------*/
/********************************* MAIN **********************************/
/*-----------------------------------------------------------------------*/


int main(int argc, char *argv[])
{

    Options opts(argc, argv);
 
    const bool getBinaryData = true;
    pwiz::msdata::FullReaderList readers;
    pwiz::msdata::MSDataFile msd(opts.mzML_file, &readers);
    pwiz::msdata::SpectrumList& spectrumList = *msd.run.spectrumListPtr;
    pwiz::msdata::SpectrumPtr spectrum;
    std::vector<pwiz::msdata::MZIntensityPair> mz_mu_pairs;
    pwiz::msdata::MZIntensityPair pair;
    
    float  rt_sigma     = opts.rt_width / 2.355;
    double mz_ppm_sigma = opts.mz_width / 2.355e6;
    int    rt_len       = spectrumList.size();
    int    mid_win      = rt_len / 2;
    pwiz::msdata::SpectrumPtr mz_mu_vect = spectrumList.spectrum(mid_win, 
                                                            getBinaryData);
    double lo_tol = 1.0 - opts.mz_sigma * mz_ppm_sigma;
    double hi_tol = 1.0 + opts.mz_sigma * mz_ppm_sigma;

    std::cout << "RT Sigma: " << rt_sigma     << std::endl
              << "MZ PPM:   " << mz_ppm_sigma << std::endl
              << "RT Len:   " << rt_len       << std::endl
              << "Mid Win:  " << mid_win      << std::endl
              << "Lo Tol:   " << lo_tol       << std::endl
              << "Hi Tol:   " << hi_tol       << std::endl;
    
    std::vector<double> points_lo_lo;
    std::vector<double> points_lo_hi;
    std::vector<double> points_hi_lo;
    std::vector<double> points_hi_hi;

    std::vector<std::vector<double>> data_lo;
    std::vector<std::vector<double>> data_hi;
    std::vector<std::vector<double>> shape_lo;
    std::vector<std::vector<double>> shape_hi;
    std::vector<int> len_lo;
    std::vector<int> len_hi;

    mz_mu_vect->getMZIntensityPairs(mz_mu_pairs);
    for (auto pair : mz_mu_pairs) {
        points_lo_lo.push_back(pair.mz * lo_tol);
        points_lo_hi.push_back(pair.mz * hi_tol);
        points_hi_lo.push_back((pair.mz + opts.mz_delta) * lo_tol);
        points_hi_hi.push_back((pair.mz + opts.mz_delta) * hi_tol);

        std::vector<double> data;
        data_lo.push_back(data);
        data_hi.push_back(data);
        shape_lo.push_back(data);
        shape_hi.push_back(data);
        len_lo.push_back(0);
        len_hi.push_back(0);
    }

    std::cout.precision(10);
    std::cout << "Points LoLo: " << points_lo_lo[0]    << std::endl
              << "Points LoHi: " << points_lo_hi[0]    << std::endl
              << "Points HiLo: " << points_hi_lo[0]    << std::endl
              << "Points HiHI: " << points_hi_hi[0]    << std::endl
              << "Data Lo:     " << data_lo[0].size()  << std::endl
              << "Data Hi:     " << data_hi[0].size()  << std::endl
              << "Shape Lo:    " << shape_lo[0].size() << std::endl
              << "Shape Hi:    " << shape_hi[0].size() << std::endl
              << "Len Lo:      " << len_lo[0]          << std::endl
              << "Len Hi:      " << len_hi[0]          << std::endl;   
    
    std::vector<float> rt_shape;

    for (int i = 0; i < rt_len; ++i) {

        float pt = (i - mid_win) / rt_sigma;
        pt = -0.5 * pt * pt;
        pt = exp(pt) / (rt_sigma * root2pi); 
    
        rt_shape.push_back(pt);
    }

    std::cout << "RT Shape: " << rt_shape.size() << std::endl
              << "RT Shape: " << rt_shape[0]     << std::endl;

    for (size_t mzi = 0; mzi < mz_mu_pairs.size(); ++mzi) {
    
        double lo_tol_lo = points_lo_lo[mzi];
        double lo_tol_hi = points_lo_hi[mzi];
        double hi_tol_lo = points_hi_lo[mzi];
        double hi_tol_hi = points_hi_hi[mzi];
        double centre    = mz_mu_pairs[mzi].mz;
        double sigma     = centre * mz_ppm_sigma;

        //std::cout << "CENTRE: " << centre << std::endl;

        pwiz::analysis::SpectrumList_MZWindow lo_window(
                                                    msd.run.spectrumListPtr,
                                                    lo_tol_lo, lo_tol_hi);
        pwiz::analysis::SpectrumList_MZWindow hi_window(
                                                    msd.run.spectrumListPtr,
                                                    hi_tol_lo, hi_tol_hi);
            
        //std::cout << "CENTRE2: " << centre << std::endl;
    
        for (int rowi = 0; rowi < rt_len; ++rowi) {
        
            //std::cout << "CENTRE3: " << centre << std::endl;
           
            centre = mz_mu_pairs[mzi].mz;
            
            pwiz::msdata::SpectrumPtr lo_spectrum;
            pwiz::msdata::SpectrumPtr hi_spectrum;
            std::vector<pwiz::msdata::MZIntensityPair> lo_pairs;
            std::vector<pwiz::msdata::MZIntensityPair> hi_pairs;
            
            lo_spectrum = lo_window.spectrum(rowi, getBinaryData);
            hi_spectrum = hi_window.spectrum(rowi, getBinaryData);
        
            lo_spectrum->getMZIntensityPairs(lo_pairs);
            hi_spectrum->getMZIntensityPairs(hi_pairs);
        
            float rt_lo = rt_shape[rowi];
            float rt_hi = rt_lo;

            //std::cout << "CENTRE4: " << centre << std::endl;
            if (lo_pairs.size() > 0) {
            
                for (auto pair : lo_pairs) {
                    //std::cout.precision(10);
                    //std::cout << "Centre: " << centre << std::endl;
                    //std::cout << "Sigma: " << sigma << std::endl;
                    //std::cout << "RT Lo: " << rt_lo << std::endl;
                    //std::cout << "MZ: " << pair.mz << std::endl;
                    float mz = (pair.mz - centre) / sigma;
                    //std::cout << "MZ: " << mz << std::endl;
                    //if ((-0.5 * mz * mz) < -183506) { 
                    //    std::cout << "MZ i: " << mzi << std::endl;
                    //    std::cout << "Lo Tol: " << lo_tol_lo << std::endl;
                    //    std::cout << "Hi Tol: " << lo_tol_hi << std::endl;
                    //    std::cout << "Row i: " << rowi << std::endl;
                    //    std::cout << "Centre: " << centre << std::endl;
                    //    std::cout << "Sigma:  " << sigma << std::endl;
                    //    std::cout << "Pair mz: " << pair.mz << std::endl;
                    //    std::cout << "MZ: " << mz << std::endl;
                    //}
                    
                    mz = -0.5 * mz * mz;
                    //std::cout << "MZ: " << mz << std::endl;
                    //if (exp(mz) == 0) {
                    //    std::cout << "ZERO EXP" << std::endl;
                    //    std::cout << mz << std::endl;
                    //    exit(1);
                    //}
                    mz = rt_lo * exp(mz) / (sigma * root2pi);
                    //if (mz == 0) {
                    //    std::cout << "ZERO POINT" << std::endl;
                    //    std::cout << rt_lo << std::endl;
                    //    std::cout << sigma << std::endl;
                    //    std::cout << root2pi << std::endl;
                    //    exit(1);
                    //}
                    //std::cout << "MZ: " << mz << std::endl;
                    shape_lo[mzi].push_back(mz);
                    data_lo[mzi].push_back(pair.intensity);
                    
                }
                len_lo[mzi] += lo_pairs.size();

            } else {
                data_lo[mzi].push_back(0);
                shape_lo[mzi].push_back(rt_lo / (sigma * root2pi));
            }
            
            centre += opts.mz_delta;
            sigma = centre * mz_ppm_sigma;
            
            if (hi_pairs.size() > 0) {
            
                for (auto pair : hi_pairs) {
                    float mz = (pair.mz - centre) / sigma;
                    mz = -0.5 * mz * mz;
                    mz = rt_hi * exp(mz) / (sigma * root2pi);
                    shape_hi[mzi].push_back(mz);
                    data_hi[mzi].push_back(pair.intensity);
                }
                len_hi[mzi] += hi_pairs.size();

            } else {
                data_hi[mzi].push_back(0);
                shape_hi[mzi].push_back(rt_hi / (sigma * root2pi));
            }
        }
    }

    std::cout << "Data Lo:     " << data_lo.size()     << std::endl
              << "Data Hi:     " << data_hi.size()     << std::endl
              << "Data Lo 0:   " << data_lo[0].size()  << std::endl
              << "Data Hi 0:   " << data_hi[0].size()  << std::endl
              << "Data Lo 0:   " << data_lo[0][0]      << std::endl
              << "Data Hi 0:   " << data_hi[0][0]      << std::endl
              << "Shape Lo:    " << shape_lo.size()    << std::endl
              << "Shape Hi:    " << shape_hi.size()    << std::endl
              << "Shape Lo 0:  " << shape_lo[0].size() << std::endl
              << "Shape Hi 0:  " << shape_hi[0].size() << std::endl
              << "Shape Lo 0:  " << shape_lo[0][0]     << std::endl
              << "Shape Hi 0:  " << shape_hi[0][0]     << std::endl
              << "Len Lo:      " << len_lo[0]          << std::endl
              << "Len Hi:      " << len_hi[0]          << std::endl;   
    
    //for (auto v : shape_lo[0]) {
    //    std::cout << v << std::endl;
    //}

    for (size_t leni = 0; leni < len_lo.size(); ++leni) {
        if (len_lo[leni] < opts.min_sample) {
            data_lo[leni]  = {0.0};
            shape_lo[leni] = {0.0};
        }
    }

    for (size_t leni = 0; leni < len_hi.size(); ++leni) {
        if (len_hi[leni] < opts.min_sample) {
            data_hi[leni]  = {0.0};
            shape_hi[leni] = {0.0};
        } else {
            for (auto& s : shape_hi[leni]) {
                s *= opts.intensity_ratio;
            }
        }
    }
    
    std::cout << "Data Lo:     " << data_lo.size()     << std::endl
              << "Data Hi:     " << data_hi.size()     << std::endl
              << "Data Lo 0:   " << data_lo[0].size()  << std::endl
              << "Data Hi 0:   " << data_hi[0].size()  << std::endl
              << "Data Lo 0:   " << data_lo[0][0]      << std::endl
              << "Data Hi 0:   " << data_hi[0][0]      << std::endl
              << "Shape Lo:    " << shape_lo.size()    << std::endl
              << "Shape Hi:    " << shape_hi.size()    << std::endl
              << "Shape Lo 0:  " << shape_lo[0].size() << std::endl
              << "Shape Hi 0:  " << shape_hi[0].size() << std::endl
              << "Shape Lo 0:  " << shape_lo[0][0]     << std::endl
              << "Shape Hi 0:  " << shape_hi[0][0]     << std::endl;


    std::vector<std::vector<double>> dataAB;
    std::vector<double> nAB;

    for (size_t i = 0; i < data_lo.size(); ++i) {
        
        std::vector<double> dataAB_row;
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

    std::cout << "Data AB:   " << dataAB.size()    << std::endl
              << "Data AB 0: " << dataAB[0].size() << std::endl
              << "Data AB 0: " << dataAB[0][0]     << std::endl
              << "nAB:       " << nAB.size()       << std::endl
              << "nAB:       " << nAB[0]           << std::endl;
     
    std::vector<std::vector<double>> shapeAB;
    std::vector<std::vector<double>> shapeA0;
    std::vector<std::vector<double>> shapeB0;
    std::vector<std::vector<double>> shape1r;

    for (size_t i = 0; i < shape_lo.size(); ++i) {
        
        std::vector<double> shapeAB_row;
        std::vector<double> shapeA0_row;
        std::vector<double> shapeB0_row;
        std::vector<double> shape1r_row;
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

    std::cout << "Shape AB:   " << shapeAB.size()    << std::endl
              << "Shape AB 0: " << shapeAB[0].size() << std::endl
              << "Shape AB 0: " << shapeAB[0][0]     << std::endl
              << "Shape A0:   " << shapeA0.size()    << std::endl
              << "Shape A0 0: " << shapeA0[0].size() << std::endl
              << "Shape A0 0: " << shapeA0[0][0]     << std::endl
              << "Shape B0:   " << shapeB0.size()    << std::endl
              << "Shape B0 0: " << shapeB0[0].size() << std::endl
              << "Shape B0 0: " << shapeB0[0][0]     << std::endl
              << "Shape 1r:   " << shape1r.size()    << std::endl
              << "Shape 1r 0: " << shape1r[0].size() << std::endl
              << "Shape 1r 0: " << shape1r[0][0]     << std::endl;

    //for (auto s : shapeAB[0]) {
    //    std::cout << s << std::endl;
    //}

    dataAB  = apply_vect_func(dataAB,  centre_vector);
    shapeAB = apply_vect_func(shapeAB, centre_vector);
    shapeA0 = apply_vect_func(shapeA0, centre_vector);
    shapeB0 = apply_vect_func(shapeB0, centre_vector);
    shape1r = apply_vect_func(shape1r, centre_vector);

    std::cout << "Data AB 0:  " << dataAB[0][0]      << std::endl
              << "Shape AB 0: " << shapeAB[0][0]     << std::endl
              << "Shape A0 0: " << shapeA0[0][0]     << std::endl
              << "Shape B0 0: " << shapeB0[0][0]     << std::endl
              << "Shape 1r 0: " << shape1r[0][0]     << std::endl;
/*               
    std::vector<std::vector<double>> data2AB;
    std::vector<std::vector<double>> shape2AB;
    std::vector<std::vector<double>> shape2A0;
    std::vector<std::vector<double>> shape2B0;
    std::vector<std::vector<double>> shape21r;

    data2AB  = apply_vect_func(dataAB,  square_vector);
    shape2AB = apply_vect_func(shapeAB, square_vector);
    shape2A0 = apply_vect_func(shapeA0, square_vector);
    shape2B0 = apply_vect_func(shapeB0, square_vector);
    shape21r = apply_vect_func(shape1r, square_vector);

    std::vector<double> SSY;
    std::vector<double> SSXAB;
    std::vector<double> SSXA0;
    std::vector<double> SSXB0;
    std::vector<double> SSX1r;

    SSY   = reduce_2D_vect(data2AB,  sum_vector);
    SSXAB = reduce_2D_vect(shape2AB, sum_vector);
    SSXA0 = reduce_2D_vect(shape2A0, sum_vector);
    SSXB0 = reduce_2D_vect(shape2B0, sum_vector);
    SSX1r = reduce_2D_vect(shape21r, sum_vector);

    std::vector<std::vector<double>> datashape;
    std::vector<double> SXYAB;
    std::vector<double> SXYA0;
    std::vector<double> SXYB0;
    std::vector<double> SXY1r;
    std::vector<double> SXYABA0;
    std::vector<double> SXYABB0;
    std::vector<double> SXYAB1r;

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

    std::vector<double> correlAB;
    std::vector<double> correlA0;
    std::vector<double> correlB0;
    std::vector<double> correl1r;
    std::vector<double> correlABA0;
    std::vector<double> correlABB0;
    std::vector<double> correlAB1r;

    correlAB   = correl_vectors(SXYAB,   SSXAB, SSY);
    correlA0   = correl_vectors(SXYA0,   SSXA0, SSY);
    correlB0   = correl_vectors(SXYB0,   SSXB0, SSY);
    correl1r   = correl_vectors(SXY1r,   SSX1r, SSY);
    correlABA0 = correl_vectors(SXYABA0, SSXAB, SSXA0);
    correlABB0 = correl_vectors(SXYABB0, SSXAB, SSXB0);
    correlAB1r = correl_vectors(SXYAB1r, SSXAB, SSX1r);
    
    std::vector<double> rm2ABA0;
    std::vector<double> rm2ABB0;
    std::vector<double> rm2AB1r;

    rm2ABA0 = rm_vectors(correlAB, correlA0);
    rm2ABB0 = rm_vectors(correlAB, correlB0);
    rm2AB1r = rm_vectors(correlAB, correl1r);
   
    std::vector<double> fABA0;
    std::vector<double> fABB0;
    std::vector<double> fAB1r;

    fABA0 = f_vectors(correlABA0, rm2ABA0);
    fABB0 = f_vectors(correlABB0, rm2ABB0);
    fAB1r = f_vectors(correlAB1r, rm2AB1r);

    std::vector<double> hABA0;
    std::vector<double> hABB0;
    std::vector<double> hAB1r;

    hABA0 = h_vectors(fABA0, rm2ABA0);
    hABB0 = h_vectors(fABB0, rm2ABB0);
    hAB1r = h_vectors(fAB1r, rm2AB1r);

    std::for_each(nAB.begin(), nAB.end(), [](double& d) { d-=3.0;});
    std::transform(nAB.begin(), nAB.end(), nAB.begin(), 
                                                 (double(*)(double)) sqrt);
   

    std::vector<double> zABA0;
    std::vector<double> zABB0;
    std::vector<double> zAB1r;

    zABA0 = z_vectors(correlAB, correlA0, nAB, correlABA0, hABA0);
    zABB0 = z_vectors(correlAB, correlB0, nAB, correlABB0, hABB0);
    zAB1r = z_vectors(correlAB, correl1r, nAB, correlAB1r, hAB1r);
    
    std::vector<double> min_score;

    for (size_t idx = 0; idx < zABA0.size(); ++idx) {
        double zA0 = zABA0[idx];
        double zB0 = zABB0[idx];
        double z1r = zAB1r[idx];
        min_score.push_back(std::min({zA0, zB0, z1r, 0.0}));
    }

    std::vector<std::vector<double>> score = {min_score, correlAB, correlA0,
                                              correlB0, correl1r};

    pwiz::msdata::SpectrumInfo spectrum_info;
    spectrum_info.update(*mz_mu_vect, getBinaryData);
    double rt = spectrum_info.retentionTime;

    std::ofstream outfile;
    outfile.open("output.txt");

    for (size_t idx = 0; idx < mz_mu_pairs.size(); ++idx) {
        double mz  = mz_mu_pairs[idx].mz;
        double amp = mz_mu_pairs[idx].intensity;
        double ms  = min_score[idx];
        double AB  = correlAB[idx];
        double A0  = correlA0[idx];
        double B0  = correlB0[idx];
        double r1  = correl1r[idx];

        outfile << rt << ", " << mz << ", " << amp << ", " << ms << ", "  
                << AB << ", " << A0 << ", " << B0  << ", " << r1 << "\n"; 
    }

    outfile.close();
*/
    std::cout << "Done!" << std::endl;
    return 0;
}


/*-----------------------------------------------------------------------*/
/************************* FUNCTION DEFINITIONS **************************/
/*-----------------------------------------------------------------------*/


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
    cout                                                            << endl;
    cout << "arguments: " << "mzML_file     path to mzML file"      << endl;
    cout                                                            << endl;
    cout << "example:   " << cmd << " example.mzML"                 << endl;
    cout                                                            << endl;
}

std::vector<double> centre_vector(std::vector<double> vect)
{
    double sum  = std::accumulate(vect.begin(), vect.end(), 0.0);
    double mean = sum / vect.size();
    //std::cout << "Sum: " << sum << std::endl;
    //std::cout << "Mean: " << mean << std::endl;
    //exit(1);
    std::vector<double> centered;

    for (auto v : vect) {
        centered.push_back(v - mean);
    }

    return centered;
}

std::vector<double> square_vector(std::vector<double> vect)
{
    std::vector<double> squared;

    for (auto v : vect) {
        squared.push_back(v * v);
    }

    return squared;
}

double sum_vector(std::vector<double> vect)
{
    double sum = std::accumulate(vect.begin(), vect.end(), 0.0);

    return sum;
}
   
std::vector<double> mult_vectors(std::vector<double> vect1, 
                                                std::vector<double> vect2)
{
    std::vector<double> mult;

    if (vect1.size() != vect2.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    for (size_t idx = 0; idx < vect1.size(); ++idx) {
        mult.push_back(vect1[idx] * vect2[idx]);
    }

    return mult;
}
   
std::vector<double> div_vectors(std::vector<double> vect1, 
                                                std::vector<double> vect2)
{
    std::vector<double> divided;

    if (vect1.size() != vect2.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    for (size_t idx = 0; idx < vect1.size(); ++idx) {
        divided.push_back(vect1[idx] / vect2[idx]);
    }

    return divided;

}

std::vector<double> correl_vectors(std::vector<double> vect1,
                    std::vector<double> vect2, std::vector<double> vect3)
{
    std::vector<double> correlated;
    std::vector<double> mult;

    mult = mult_vectors(vect2, vect3);
    std::transform(mult.begin(), mult.end(), mult.begin(), 
                                            (double(*)(double)) std::sqrt);
    correlated = div_vectors(vect1, mult);

    for (auto& c : correlated) {
        if(c < 0) {
            c = 0;
        }
    }

    return correlated;
}

std::vector<double> rm_vectors(std::vector<double> vect1, 
                                                std::vector<double> vect2)
{
    std::vector<double> rm;

    if (vect1.size() != vect2.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    for (size_t idx = 0; idx < vect1.size(); ++idx) {
        double tmp = (vect1[idx] * vect1[idx]) + (vect2[idx] * vect2[idx]);
        rm.push_back(0.5 * tmp);
    }

    return rm;

}

std::vector<double> f_vectors(std::vector<double> correl_vect,
                                              std::vector<double> rm_vect)
{
    std::vector<double> f_vect;

    if (correl_vect.size() != rm_vect.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    for (size_t idx = 0; idx < correl_vect.size(); ++idx) {
        double correl = correl_vect[idx];
        double rm     = rm_vect[idx];

        f_vect.push_back((1.0 - correl) / (2.0 * (1.0 - rm)));
    }

    for (auto& f : f_vect) {
        if(f > 1.0) {
            f = 1.0;
        }
    }

    return f_vect;
}

std::vector<double> h_vectors(std::vector<double> f_vect,
                                              std::vector<double> rm_vect)
{
    std::vector<double> h_vect;

    if (f_vect.size() != rm_vect.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    for (size_t idx = 0; idx < f_vect.size(); ++idx) {
        double f  = f_vect[idx];
        double rm = rm_vect[idx];

        h_vect.push_back((1.0 - f * rm) / (1.0 - rm));
    }

    return h_vect;
}

std::vector<double> z_vectors(std::vector<double> cor1,
                              std::vector<double> cor2,
                              std::vector<double> sqrtn,
                              std::vector<double> cross_cor,
                              std::vector<double> h_vect)
{
    std::vector<double> z_vect;

    for (size_t idx = 0; idx < cor1.size(); ++idx) {
        
        double z1  = std::atanh(cor1[idx]);
        double z2  = std::atanh(cor2[idx]);
        
        double num   = (z1 - z2) * sqrtn[idx];
        double denom = 2.0 * (1.0 - cross_cor[idx]) * h_vect[idx];

        z_vect.push_back(num / std::sqrt(denom));
    }

    return z_vect;
}

template <typename T, typename F>
std::vector<T> apply_vect_func(std::vector<T> vect, F func)
{
    std::vector<T> applied;
    
    for (auto v : vect) {
        applied.push_back(func(v));
    }

    return applied;
}

template <typename T, typename F>
std::vector<T> apply_vect_func(std::vector<T> vect1, std::vector<T> vect2, 
                                                                     F func)
{
    std::vector<T> applied;
    
    if (vect1.size() != vect2.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    for (size_t idx = 0; idx < vect1.size(); ++idx) {
        applied.push_back(func(vect1[idx], vect2[idx]));
    }

    return applied;
}

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
    mzML_file       = "";

    // Show usage and exit if no options are given
    if (argc == 1) {
        show_usage(argv[0]);
        exit(1);
    }

    while ((opt = getopt(argc, argv, "hd:i:r:R:p:m:M:D:s:")) != -1){
        
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
        }
    }

    for (opt_idx = optind; opt_idx < argc; opt_idx++) {

        if (mzML_file == "") { 
            mzML_file = argv[opt_idx];
        } else {
            std::cout << "Too many arguments supplied. See usage.";
            std::cout << std::endl;
            exit(1);
        }
    }

    if (mzML_file == "") {
        std::cout << "Insufficient arguments supplies. See usage.";
        std::cout << std::endl;
        exit(1);
    }
}


/*-----------------------------------------------------------------------*/
/******************************* OLD CODE ********************************/
/*-----------------------------------------------------------------------*/


