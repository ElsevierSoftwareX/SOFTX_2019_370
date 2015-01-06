#include <iostream>
#include <stdlib.h>

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

    cout << "Usage:   " << cmd <<" [-option] [argument]"<<endl;
    cout << "option:  " << "-h  show help information"<<endl;
    cout << "         " << "-u username"<<endl;
    cout << "         " << "-p  password"<<endl;
    cout << "         " << "-s  save the password: 0(save password) 1(forget password)"<<endl;
    cout << "         " << "-v  show version infomation"<<endl;
    cout << "example: " << cmd <<" -uusername -ppassword -s1"<<endl;
}

int main(int argc, char *argv[])
{
 
    if (argc == 1) {
        show_usage(argv[0]);
        exit(1);
    } else { 
        /*for loop,each loop print a argument once a time. Note that the 
          loop begin with argv[1], because argv[0] represent the program's 
          name.*/
        std::cout<<"The arguments you put are:"<<std::endl;
        for(int i=1;i<argc;i++) {
            std::cout<<argv[i]<<std::endl;
        }
        return 0;
    }

/*
    std::string mz_file = "./data/testing_mz.npy";
    arma::mat mz_mat = load_npy_file(mz_file);
    mz_mat.submat(0, 0, 5, 5).print("MZ Data:");

    std::string time_file = "./data/testing_time.npy";
    arma::mat time_mat = load_npy_file(time_file);
    time_mat.submat(0, 0, 5, 5).print("Time Data:");

    std::string intensity_file = "./data/testing_intensity.npy";
    arma::mat intensity_mat = load_npy_file(intensity_file);
    intensity_mat.submat(0, 0, 5, 5).print("Intensity Data:");
*/
    
    return 0;
}


