#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>
#include "constants.h"
#include "options.h"

namespace po = boost::program_options;
using namespace std;

Options::Options(int argc, char* argv[])
{
    intensity_ratio = default_intensity_ratio;
    rt_width = default_rt_width;
    rt_sigma = default_rt_sigma;
    ppm = default_ppm;
    mz_width = default_fwhm;
    mz_sigma = default_mz_sigma;
    mz_delta = default_mz_delta;
    min_sample = default_min_sample;
    in_file = "";
    out_file = "";
    debug = false;
    num_threads = 1;

    string iratio_str = "Ratio of doublet intensities (isotope / parent). Defaults to " + to_string(default_intensity_ratio);
    string rtwidth_str = "Full width at half maximum for retention time in number of scans. Defaults to " + to_string(default_rt_width);
    string rtwindow_str = "Retention time width boundary in standard deviations. Defaults to " + to_string(default_rt_sigma);
    string ppm_str = "M/Z tolerance in parts per million. Defaults to " + to_string(default_ppm);
    string mzwidth_str = "M/Z full width at half maximum in parts per million. Defaults to " + to_string(default_fwhm);
    string mzsigma_str = "M/Z window boundary in standard deviations. Defaults to " + to_string(default_mz_sigma);
    string mzdelta_str = "M/Z delta for doublets. Defaults to " + to_string(default_mz_delta);
    string minsample_str = "Minimum number of data points required in each sample region. Defaults to " + to_string(default_min_sample);
    string threads_str = "Number of threads to use. Defaults to "  + to_string(num_threads);

    po::options_description desc(program_name + " allowed options");
    desc.add_options()
        ("help,h", "Show this help information")
        ("iratio,a", po::value<double>(), iratio_str.c_str())
        ("rtwidth,r", po::value<double>(), rtwidth_str.c_str())
        ("rtwindow,t", po::value<double>(), rtwindow_str.c_str())
        ("ppm,p", po::value<double>(), ppm_str.c_str())
        ("mzwidth,m", po::value<double>(), mzwidth_str.c_str())
        ("mzwindow,z", po::value<double>(), mzsigma_str.c_str())
        ("mzdelta,d", po::value<double>(), mzdelta_str.c_str())
        ("mindata,n", po::value<int>(), minsample_str.c_str())
        ("debug", "Generate debugging output")
        ("threads,j", po::value<int>(), "Number of threads to use")
        ("infile,i", po::value<string>()->required(), "Input mzML file")
        ("outfile,o", po::value<string>()->required(), "Output mzML file");

    if (argc <= 1)
    {
        cout << program_name << " insufficient command line arguments." << endl;
        cout << desc << "\n";
        exit(0);
    }

    po::variables_map vm;

    try
    {
       po::store(po::parse_command_line(argc, argv, desc), vm);

       if (vm.count("help"))
       {
          cout << desc << "\n";
          exit(0);
       }
       po::notify(vm); // throws on error, so do after help in case
                       // there are any problems
     }
     catch(boost::program_options::required_option& e)
     {
          std::cerr << program_name << " ERROR: " << e.what() << std::endl << std::endl;
          exit(-1);
     }
     catch(boost::program_options::error& e)
     {
          std::cerr << program_name << " ERROR: " << e.what() << std::endl << std::endl;
          exit(-1);
     }

    if (vm.count("help")) {
        cout << desc << "\n";
        exit(0);
    }
    if (vm.count("iratio")) {
        intensity_ratio = vm["iratio"].as<double>();
    }
    if (vm.count("rtwidth")) {
        rt_width = vm["rtwidth"].as<double>();
    }
    if (vm.count("rtwindow")) {
        rt_sigma = vm["rtwindow"].as<double>();
    }
    if (vm.count("ppm")) {
        ppm = vm["ppm"].as<double>();
    }
    if (vm.count("mzwidth")) {
        mz_width = vm["mzwidth"].as<double>();
    }
    if (vm.count("mzwindow")) {
        mz_sigma = vm["mzwindow"].as<double>();
    }
    if (vm.count("mzdelta")) {
        mz_delta = vm["mzdelta"].as<double>();
    }
    if (vm.count("mindata")) {
        min_sample = vm["mindata"].as<double>();
    }
    if (vm.count("infile")) {
        in_file = vm["infile"].as<string>();
    }
    if (vm.count("outfile")) {
        out_file = vm["outfile"].as<string>();
    }
    if (vm.count("debug")) {
        debug = true;
    }
    if (vm.count("threads")) {
        int requested_threads = vm["threads"].as<int>();
        if (requested_threads < 1)
        {
            cerr << program_name << " ERROR: number of requested threads may not be less than 1";
            exit(-1);
        }
        num_threads = requested_threads;
    }
}
