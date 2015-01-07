#include <iostream>
#include <unistd.h>

#include <armadillo>
#include "Numpy.hpp"

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

/**
  Read in a NumPy array and return an Armadillo matrix
*/
arma::mat load_npy_file(const std::string& filename)
{
    std::vector<int> shape;
    std::vector<double> data;
    
    aoba::LoadArrayFromNumpy(filename, shape, data);

    int nrows = shape[0];
    int ncols = shape[1];
    int row = 0;
    int col = 0;
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

void show_usage(char *cmd)
{
    using namespace std;

    cout << "Usage:     " << cmd << " [-options] [arguments]"       << endl;
    cout                                                            << endl;
    cout << "options:   " << "-h  show this help information"       << endl;
    cout << "           " << "-d  directory  path to data files"    << endl;
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
    cout << "arguments: " << "mz_file    name of mz NumPy file"     << endl;
    cout << "           " << "time_file  name of time NumPy file"   << endl;
    cout << "           " << "int_file   name of intensity NumPy file"     ;
    cout                                                            << endl;
    cout                                                            << endl;
    cout << "example:   " << cmd << " -d ./data/ mz.npy time.npy "         ;
    cout << "intensity.npy"                                         << endl;
    cout                                                            << endl;
}

int main(int argc, char *argv[])
{
 
    char opt;
    int opt_idx;

    float intensity_ratio = default_intensity_ratio;
    float rt_width        = default_rt_width;
    float rt_sigma        = default_rt_sigma;
    float ppm             = default_ppm;
    float mz_width        = default_fwhm;
    float mz_sigma        = default_mz_sigma;
    float mz_delta        = default_mz_delta;
    float min_sample      = default_min_sample;

    std::string file_dir = "./";
    std::string mz_file = "";
    std::string time_file = "";
    std::string intensity_file = "";

    // Show usage if no options are given
    if (argc == 1) {
        show_usage(argv[0]);
        exit(1);
    }

    while ((opt = getopt (argc, argv, "hd:i:r:R:p:m:M:D:s:")) != -1) {
        
        switch (opt) {
            case 'h':
                show_usage(argv[0]);
                exit(1);
                break;
            case 'd':
                file_dir = optarg;
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
        
        if (mz_file == "") { 
            mz_file = file_dir + argv[opt_idx];
        } else if (time_file == "") {
            time_file = file_dir + argv[opt_idx];
        } else if (intensity_file == "") {
            intensity_file = file_dir + argv[opt_idx];
        } else {
            std::cout << "Too many arguments supplied. See usage.";
            std::cout << std::endl;
            exit(1);
        }

    }

    if (intensity_file == "") {
        std::cout << "Insufficient arguments supplies. See usage.";
        std::cout << std::endl;
        exit(1);
    }
    
    std::cout << "File Directory: " << file_dir << std::endl;
    std::cout << "MZ File: " << mz_file << std::endl;
    std::cout << "Time File: " << time_file << std::endl;
    std::cout << "Intensity File: " << intensity_file << std::endl;
    
    arma::mat mz_mat = load_npy_file(mz_file);
    mz_mat.submat(0, 0, 5, 5).print("MZ Data:");

    arma::mat time_mat = load_npy_file(time_file);
    time_mat.submat(0, 0, 5, 5).print("Time Data:");

    arma::mat intensity_mat = load_npy_file(intensity_file);
    intensity_mat.submat(0, 0, 5, 5).print("Intensity Data:");
    
    return 0;
}


