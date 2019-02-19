#include <iostream>
#include <iterator>
#include "constants.h"
#include "options.h"
#include "cxxopts.h"
#include "version.h"

using namespace std;

Options::Options(int argc, char* argv[])
{
    intensity_ratio = default_intensity_ratio;
    mz_lower = 0;
    mz_upper = 0;  // use 0 to show not set
    confidence = 0;
    in_file = "";
    out_file = "";
    debug = false;
    num_threads = 1;
    input_spectrum_cache_size = default_input_spectrum_cache_size;
    int num_args;

    string mzlower_str = "Lower M/Z offset for local max window, e.g. 0.01. Must be used with '--mzupper', can not be used with '--mzwidth'";
    string mzupper_str = "Upper M/Z offset for local max window, e.g. 0.5. Must be used with '--mzlower', can not be used with '--mzwidth'";
    string iratio_str = "Ratio of doublet intensities (isotope / parent). Defaults to " + to_string(default_intensity_ratio);
    string rtwidth_str = "Full width at half maximum for retention time in number of scans. Eg: 10";
    string mzwidth_str = "M/Z full width at half maximum in parts per million. Eg: 230. Must be used with '--mzdelta', can not be used with '--mzlower' and '--mzupper'.";
    string mzdelta_str = "M/Z delta for doublets. Eg: 6.0201. Must be used with '--rtwidth' and '--mzwidth', can not be used with '--mzlower' and '--mzupper'.";
    string confidence_str = "Lower confidence interval to apply during scoring (In standard deviations, e.g. 1.96 for a 95% CI). Default: ignore confidence intervals";
    string threads_str = "Number of threads to use. Defaults to "  + to_string(num_threads);
    string desc = "Detect twin ion signal in Mass Spectrometry data.\n"
                  "With 'rtwidth', 'mzwidth', and 'mzdelta', score data based on how well it matches twin ion signal\n"
                  "With 'rtwidth', 'mzlower', and 'mzupper', output list of local maxima in bounds.";
    string input_spectrum_cache_size_str = "Number of input spectra to retain in cache. Defaults to " + to_string(default_input_spectrum_cache_size);

    string msgs = "";
    try {
        cxxopts::Options options("HiTIME-CPP", desc);
        options.add_options()
            ("h,help", "Show this help information.")
            ("l,mzlower", mzlower_str, cxxopts::value<double>())
            ("u,mzupper", mzupper_str, cxxopts::value<double>())
            ("a,iratio", iratio_str, cxxopts::value<double>())
            ("r,rtwidth", rtwidth_str, cxxopts::value<double>())
            ("m,mzwidth", mzwidth_str, cxxopts::value<double>())
            ("d,mzdelta", mzdelta_str, cxxopts::value<double>())
            ("z,confidence", confidence_str, cxxopts::value<double>())
            ("debug", "Generate debugging output")
            ("version", "Print version number and exit")
            ("j,threads", threads_str, cxxopts::value<int>())
            ("c,cache", input_spectrum_cache_size_str , cxxopts::value<int>())
            ("i,infile", "Input mzML file", cxxopts::value<string>())
            ("o,outfile", "Output mzML file", cxxopts::value<string>());

        num_args = argc;
        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            cout << options.help() << endl;
            exit(0);
        }
        if (result.count("version")) {
            cout << " version " << HITIME_VERSION << endl; 
            exit(0);
        }

        if (num_args <= 1)
        {
            msgs += " insufficient command line arguments.\n";
        }

        if (result.count("mzlower")) {
            mz_lower = result["mzlower"].as<double>();
            if (mz_lower <= 0)
            {
                msgs += " ERROR: Lower M/Z bound must be greater than zero\n";
            }
            if (result.count("mzwidth") or result.count("mzdelta")) {
                msgs += " ERROR: 'mzlower' used with incompatable options\n";
            }
        }
        if (result.count("mzupper")) {
            mz_upper = result["mzupper"].as<double>();
            if (mz_upper <= 0)
            {
                msgs += " ERROR: Upper M/Z bound must be greater than zero\n";
            }
            if (result.count("mzwidth") or result.count("mzdelta")) {
                msgs += " ERROR: 'mzupper' used with incompatable options\n";
            }
        }
        if (result.count("iratio")) {
            intensity_ratio = result["iratio"].as<double>();
            if (intensity_ratio <= 0)
            {
                msgs += " ERROR: Intensity ration must be greater than zero\n";
            }
        }
        if (result.count("rtwidth")) {
            rt_width = result["rtwidth"].as<double>();
            if (rt_width <= 0)
            {
                msgs += " ERROR: Retention time full width at half maximum must be greater than zero\n";
            }
        }
        else
        {
                msgs += " ERROR: Retention time full width at half maximum is required\n";
        }
        if (result.count("mzwidth")) {
            mz_width = result["mzwidth"].as<double>();
            if (mz_width <= 0)
            {
                msgs += " ERROR: m/z full width at half maximum must be greater than zero\n";
            }
            if (result.count("mzlower") or result.count("mzupper")) {
                msgs += " ERROR: 'mzwidth' used with incompatable options\n";
            }
        }
        if (result.count("mzdelta")) {
            mz_delta = result["mzdelta"].as<double>();
            if (mz_delta <= 0)
            {
                msgs += " ERROR: m/z twin ion mass difference must be greater than zero\n";
            }
            if (result.count("mzlower") or result.count("mzupper")) {
                msgs += " ERROR: 'mzdelta' used with incompatable options\n";
            }
        }
        if (result.count("confidence")) {
            confidence = result["confidence"].as<double>();
            if (confidence <= 0)
            {
                msgs += " ERROR: confidance must be greater than zero\n";
            }
        }
        if (result.count("infile")) {
            in_file = result["infile"].as<string>();
        }
        else
        {
            msgs += "ERROR: input file name is required\n";
        }
        if (result.count("outfile")) {
            out_file = result["outfile"].as<string>();
        }
        else
        {
            msgs += "ERROR: output file name is required\n";
        }
        if (result.count("debug")) {
            debug = true;
        }
        if (result.count("threads")) {
            int requested_threads = result["threads"].as<int>();
            if (requested_threads < 1)
            {
                msgs += " ERROR: number of requested threads may not be less than 1\n";
            }
            num_threads = requested_threads;
        }
        if (result.count("cache")) {
            int requested_size = result["cache"].as<int>();
            if (requested_size < 0)
            {
                msgs += " ERROR: requested cache size must be non-negative\n";
            }
            input_spectrum_cache_size = requested_size;
        }
        if (not ((result.count("mzdelta") and result.count("mzwidth")) or
                  result.count("mzlower") and result.count("mzupper"))) {
            msgs += " ERROR: 'mzdelta' and 'mzwidth' OR 'mzlower' and 'mzupper' are required\n";
        }
        if (msgs != "") {
            cerr << program_name << endl;
            cerr << msgs << endl;
            cerr << options.help() << endl;
            exit(1);
        };
    }

    catch (const cxxopts::OptionException& e)
    {
        std::cerr << program_name << std::endl;
        std::cerr << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }
}
