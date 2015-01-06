#include <iostream>
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

int main()
{
    std::string mz_file = "./data/testing_mz.npy";
    arma::mat mz_mat = load_npy_file(mz_file);
    mz_mat.submat(0, 0, 5, 5).print("MZ Data:");

    std::string time_file = "./data/testing_time.npy";
    arma::mat time_mat = load_npy_file(time_file);
    time_mat.submat(0, 0, 5, 5).print("Time Data:");

    std::string intensity_file = "./data/testing_intensity.npy";
    arma::mat intensity_mat = load_npy_file(intensity_file);
    intensity_mat.submat(0, 0, 5, 5).print("Intensity Data:");

    return 0;
}
