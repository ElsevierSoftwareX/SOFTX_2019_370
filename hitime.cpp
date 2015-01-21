#include <iostream>
#include <unistd.h>
#include <math.h>
#include <iomanip>
#include <cstdlib>
#include <string>

#include "pwiz_tools/common/FullReaderList.hpp"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/analysis/spectrum_processing/SpectrumList_MZWindow.hpp"
//#include <armadillo>
//#include "Numpy.hpp"


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


/*-----------------------------------------------------------------------*/ 
/************************* FUNCTION DECLARATIONS *************************/
/*-----------------------------------------------------------------------*/


void show_usage(char *cmd);


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
   
    pwiz::msdata::FullReaderList readers;
    pwiz::msdata::MSDataFile msd(opts.mzML_file, &readers);

    std::cout << "Timestamp: " << msd.run.startTimeStamp << std::endl;

    pwiz::msdata::SpectrumList& spectrumList = *msd.run.spectrumListPtr;
    const bool getBinaryData = true;
    size_t numSpectra = spectrumList.size();

    std::cout << "Num Spectra: " << numSpectra << std::endl;
    
    pwiz::msdata::SpectrumPtr spectrum;
    std::vector<pwiz::msdata::MZIntensityPair> pairs;
    spectrum = spectrumList.spectrum(0, getBinaryData);
    spectrum->getMZIntensityPairs(pairs);

    std::cout << "Num Pairs: " << pairs.size() << std::endl;

    pwiz::msdata::MZIntensityPair pair;
    pair = pairs[0];

    std::cout << "MZ/Int Pair: " << pair.mz << " " << pair.intensity; 
    std::cout << std::endl;

    pwiz::analysis::SpectrumList_MZWindow mz_window(msd.run.spectrumListPtr,
                                                    150.2, 150.3);
    spectrum = mz_window.spectrum(0, getBinaryData);
    spectrum->getMZIntensityPairs(pairs);
    pair = pairs[0];

    std::cout << "WINDOW 1 (150.2 - 150.3): " << std::endl;
    std::cout << "Num Pairs: " << pairs.size() << std::endl;
    std::cout << "MZ/Int Pair: " << pair.mz << " " << pair.intensity; 
    std::cout << std::endl;

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
    mzML_file        = "";

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

/*
class MZWindow {
    
    public:
        unsigned int centre_col;
        double centre_mz;
        float tolerance;
        unsigned int col_half_length;
        std::vector<int> row_lower_bound;
        std::vector<int> row_upper_bound;

    void SetBounds (arma::mat matrix, unsigned int col,
                    unsigned int half, double mz, float tol) 
    {
        centre_col      = col;
        col_half_length = half;
        centre_mz       = mz;
        tolerance       = tol;
        
        unsigned int start_col = std::max(0U, centre_col - col_half_length);
        unsigned int end_col = std::min(centre_col + col_half_length, 
                                        matrix.n_cols-1);
        arma::uvec indexes;

        for (int col_idx = start_col; col_idx <= end_col; col_idx++) {
            
            indexes = arma::find(matrix.col(col_idx) > centre_mz - tolerance
                            && matrix.col(col_idx) < centre_mz + tolerance);
            
            if (indexes.n_elem > 0) {
                row_lower_bound.push_back(indexes(0));
                row_upper_bound.push_back(indexes(indexes.n_elem-1));
            } else {
                row_lower_bound.push_back(-1);
                row_upper_bound.push_back(-1);
            }
        }
    }

    void PrintBounds ()
    {
        for (std::vector<int>::const_iterator i = row_lower_bound.begin();
             i != row_lower_bound.end(); ++i) {
            std::cout << (*i) << " ";
        }

        std::cout << std::endl;
        for (std::vector<int>::const_iterator i = row_upper_bound.begin(); 
             i != row_upper_bound.end(); ++i) {
            std::cout << (*i) << " ";
        }
        std::cout << std::endl;
    }

};

class DoubleWindow {

    public:
        MZWindow lo_window;
        MZWindow hi_window;

    void SetWindows (arma::mat matrix, unsigned int col, unsigned int half,
                     double mz, float delta, float lo_tol, float hi_tol)
    {
        lo_window.SetBounds(matrix, col, half, mz, lo_tol);
        hi_window.SetBounds(matrix, col, half, mz+delta, hi_tol);
    }

};

*/

/**
  Read in a NumPy array and return an Armadillo matrix
*/
/*
arma::mat load_npy_file(const std::string& filename)
{
    std::vector<int> shape;
    std::vector<double> data;
    
    aoba::LoadArrayFromNumpy(filename, shape, data);

    int nrows = shape[0];
    int ncols = shape[1];
    int row   = 0;
    int col   = 0;
    std::vector<double>::const_iterator i;

    arma::mat matrix(nrows, ncols);

    for (i = data.begin(); i != data.end(); i++) {
        matrix(row, col) = (*i);
        col++;
        if (col == ncols) {
            col = 0;
            row++;
        }
    }

    return matrix;
}
*/
/*
    std::cout << "File Directory: " << file_dir << std::endl;
    std::cout << "MZ File: " << mz_file << std::endl;
    std::cout << "Time File: " << time_file << std::endl;
    std::cout << "Intensity File: " << intensity_file << std::endl;
 
*/
/*
    arma::mat mz_mat = load_npy_file(mz_file);
//    mz_mat.submat(0, 0, 5, 5).print("MZ Data:");

    arma::mat time_mat = load_npy_file(time_file);
//    time_mat.submat(0, 0, 5, 5).print("Time Data:");

    arma::mat intensity_mat = load_npy_file(intensity_file);
//    intensity_mat.submat(0, 0, 5, 5).print("Intensity Data:");


//    half_rt_window = static_cast<int>(ceil(rt_sigma * rt_width / 2.355));


    arma::mat C;
    C << 0.11 << 0.01 << 0.09 << 0.10 << 0.12 << 0.08 << 0.09 << arma::endr
      << 0.19 << 0.14 << 0.23 << 0.19 << 0.28 << 0.16 << 0.24 << arma::endr
      << 0.31 << 0.26 << 0.27 << 0.34 << 0.32 << 0.28 << 0.32 << arma::endr
      << 0.43 << 0.38 << 0.41 << 0.48 << 0.38 << 0.43 << 0.41 << arma::endr
      << 0.49 << 0.54 << 0.46 << 0.53 << 0.51 << 0.50 << 0.49 << arma::endr
      << 0.62 << 0.64 << 0.58 << 0.61 << 0.66 << 0.58 << 0.63 << arma::endr
      << 0.70 << 0.68 << 0.71 << 0.70 << 0.72 << 0.69 << 0.66 << arma::endr
      << 0.76 << 0.81 << 0.79 << 0.77 << 0.83 << 0.84 << 0.81 << arma::endr
      << 0.89 << 0.92 << 0.96 << 0.89 << 0.91 << 0.95 << 0.92 << arma::endr;

    C.print("C: ");

    MZWindow window;
    //window.SetBounds(C, 2, 3, 0.15, 2);

    //arma::mat D = mz_mat.submat(20, 20, 30, 24);
    //D.print("D:");

    //MZWindow wind2;
    //wind2.SetBounds(D, 4, 3, D(3,4), 0.02);
    //wind2.PrintBounds();

    //time_mat.submat(20, 20, 30, 24).print("Time: ");

    DoubleWindow dwind;
    dwind.SetWindows(C, 3, 2, C(2, 3), 0.3, 0.11, 0.16); 
    dwind.lo_window.PrintBounds();
    dwind.hi_window.PrintBounds();

*/
