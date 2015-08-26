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

/*
#include "pwiz_tools/common/FullReaderList.hpp"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/analysis/spectrum_processing/SpectrumList_MZWindow.hpp"
#include "pwiz/data/msdata/SpectrumInfo.hpp"
#include "pwiz/data/common/cv.hpp"
*/

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
const float default_min_sample = default_rt_width * default_rt_sigma / 2.355;
const double root2pi = sqrt(2.0 * M_PI);

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


//! @brief Print program usage information.
void show_usage(char *cmd);

int main(int argc, char** argv)
{
   // Read user options
   Options opts(argc, argv);

   IndexedMzMLFileLoader imzml;
   //MzMLFile mzml;

   // load data from an indexed MzML file
   OnDiscMSExperiment<> map;
   //MSExperiment<> map;

   imzml.load(opts.mzML_file, map);
   //mzml.load(opts.mzML_file, map);

   // Get the first spectrum in memory, do some constant (non-changing) data processing
   MSSpectrum<> s = map.getSpectrum(0);
   std::cout << "There are " << map.getNrSpectra() << " spectra in the input file." << std::endl;
   std::cout << "The first spectrum has " << s.size() << " peaks." << std::endl;

   // store the (unmodified) data in a different file
   //imzml.store("Tutorial_FileIO_output.mzML", map);

   // Calculate number of spectra for each window
   // XXX fix magic number
   int half_window = ceil(opts.rt_sigma * opts.rt_width / 2.355);


/*
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
*/

   std::cout << "Done!" << std::endl;
   return 0;

}

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
/*
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
*/

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
        std::cout << "Insufficient arguments supplied. See usage.\n";
	show_usage(argv[0]);
        exit(1);
    }
}
