#include <iostream>
#include <iterator>
#include "constants.h"
#include "options.h"
#include "cxxopts.h"

using namespace std;

Options::Options(int argc, char* argv[])
{
    list_max = false;
    intensity_ratio = default_intensity_ratio;
    confidence = 0;
    in_file = "";
    out_file = "";
    debug = false;
    num_threads = 1;
    input_spectrum_cache_size = default_input_spectrum_cache_size;
    int num_args;

    string list_max_str = "Flag, only output list of local maximum in window defined by M/Z width and retention time width. Default: not set";
    string iratio_str = "Ratio of doublet intensities (isotope / parent). Defaults to " + to_string(default_intensity_ratio);
    string rtwidth_str = "REQUIRED: Full width at half maximum for retention time in number of scans. Eg: 17";
    string mzwidth_str = "REQUIRED: M/Z full width at half maximum in parts per million. Eg: 150. If '--listmax', then upper and lower M/Z offset, e.g. 0.25";
    string mzdelta_str = "REQUIRED: M/Z delta for doublets. Eg: " + to_string(default_mz_delta);
    string confidence_str = "Lower confidence interval to apply during scoring (In standard deviations, e.g. 1.96 for a 95% CI). Default: ignore confidence intervals";
    string threads_str = "Number of threads to use. Defaults to "  + to_string(num_threads);
    string desc = "Detect twin ion signal in Mass Spectrometry data";
    string input_spectrum_cache_size_str = "Number of input spectra to retain in cache. Defaults to " + to_string(default_input_spectrum_cache_size);

    string msgs = "";
    try {
        cxxopts::Options options("HiTIME-CPP", desc);
        options.add_options()
            ("h,help", "Show this help information.")
            ("l,listmax", list_max_str, cxxopts::value<bool>())
            ("a,iratio", iratio_str, cxxopts::value<double>())
            ("r,rtwidth", rtwidth_str, cxxopts::value<double>())
            ("m,mzwidth", mzwidth_str, cxxopts::value<double>())
            ("d,mzdelta", mzdelta_str, cxxopts::value<double>())
            ("z,confidence", confidence_str, cxxopts::value<double>())
            ("debug", "Generate debugging output")
            ("j,threads", threads_str, cxxopts::value<int>())
            ("c,cache", input_spectrum_cache_size_str , cxxopts::value<int>())
            ("i,infile", "Input mzML file", cxxopts::value<string>())
            ("o,outfile", "Output mzML file", cxxopts::value<string>());

        num_args = argc;
        auto result = options.parse(argc, argv);

        if (num_args <= 1)
        {
            cout << program_name << " insufficient command line arguments." << endl;
            cout << options.help() << endl;
            exit(0);
        }

        if (result.count("help")) {
            cout << options.help() << endl;
            exit(0);
        }

        if (result.count("listmax")) {
            list_max = result["listmax"].as<bool>();
        }
        if (result.count("iratio")) {
            intensity_ratio = result["iratio"].as<double>();
        }
        if (result.count("rtwidth")) {
            rt_width = result["rtwidth"].as<double>();
            if (rt_width <= 0)
            {
                cerr << program_name << " ERROR: Retention time full width at half maximum must be greater than zero";
                exit(-1);
            }
        }
        else
        {
            msgs += "MISSING: Retention time full width at half maximum must be given.\n";
        }
        if (result.count("mzwidth")) {
            mz_width = result["mzwidth"].as<double>();
            if (mz_width <= 0)
            {
                cerr << program_name << " ERROR: m/z full width at half maximum must be greater than zero";
                exit(-1);
            }
        }
        else
        {
            msgs += "MISSING: m/z full width at half maximum must be given.\n";
        }
        if (result.count("mzdelta")) {
            mz_delta = result["mzdelta"].as<double>();
            if (mz_delta <= 0)
            {
                cerr << program_name << " ERROR: m/z twin ion mass difference must be greater than zero";
                exit(-1);
            }
        }
        // Need M/Z mass difference unless listing local maxima
        else if (result.count("listmax") == false)
        {
            msgs += "MISSING: m/z twin ion mass difference must be given.\n";
        }
        if (result.count("confidence")) {
            confidence = result["confidence"].as<double>();
            if (confidence <= 0)
            {
                cerr << program_name << " ERROR: confidance must be greater than zero";
                exit(-1);
            }
        }
        if (result.count("infile")) {
            in_file = result["infile"].as<string>();
        }
        if (result.count("outfile")) {
            out_file = result["outfile"].as<string>();
        }
        if (result.count("debug")) {
            debug = true;
        }
        if (result.count("threads")) {
            int requested_threads = result["threads"].as<int>();
            if (requested_threads < 1)
            {
                cerr << program_name << " ERROR: number of requested threads may not be less than 1";
                exit(-1);
            }
            num_threads = requested_threads;
        }
        if (result.count("cache")) {
            int requested_size = result["cache"].as<int>();
            if (requested_size < 0)
            {
                cerr << program_name << " ERROR: requested cache size must be non-negative";
                exit(-1);
            }
            input_spectrum_cache_size = requested_size;
        }
        if (msgs != "") {
            cout << program_name << endl;
            cout << msgs << endl;
            cout << options.help() << endl;
            exit(-1);
        }
    }

    catch (const cxxopts::OptionException& e)
    {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }
}
