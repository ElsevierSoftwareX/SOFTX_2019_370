#include <iostream>
#include <unistd.h>

#include <armadillo>
#include "Numpy.hpp"


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

    cout << "Usage:     " << cmd << " [-options] [arguments]" << endl;
    cout << "options:   " << "-h  show this help information" << endl;
    cout << "           " << "-d  directory  path to data files" << endl;
    cout << "arguments: " << "mz_file    name of mz NumPy file" << endl;
    cout << "           " << "time_file  name of time NumPy file" << endl;
    cout << "           " << "int_file   name of intensity NumPy file"; 
    cout << endl;
    cout << "example:   " << cmd << " -d ./data/ mz.npy time.npy ";
    cout << "intensity.npy" << endl;
}

int main(int argc, char *argv[])
{
 
    char opt;
    int opt_idx;

    std::string file_dir = "./";
    std::string mz_file = "";
    std::string time_file = "";
    std::string intensity_file = "";

    // Show usage if no options are given
    if (argc == 1) {
        show_usage(argv[0]);
        exit(1);
    }

    while ((opt = getopt (argc, argv, "hd:")) != -1) {
        
        switch (opt) {
            case 'h':
                show_usage(argv[0]);
                exit(1);
                break;
            case 'd':
                file_dir = optarg;
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


