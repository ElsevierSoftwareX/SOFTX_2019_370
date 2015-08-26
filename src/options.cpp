#include <iostream>
#include "options.h"
#include "defaults.h"

using namespace std;

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
