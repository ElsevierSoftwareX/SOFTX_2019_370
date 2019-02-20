#ifndef HITIME_OPTIONS_H
#define HITIME_OPTIONS_H

#include <unistd.h>
#include <string>

class Options {

    public:
        bool debug;
        bool version;
        double intensity_ratio; //!< Intensity ratio between lo and hi peaks.
        double rt_width; //!< Retention time FWHM in scans.
        double rt_sigma; //!< Boundary for RT width in SDs.
        double mz_width; //!< MZ FWHM in PPM.
        double mz_sigma; //!< Boundary for MZ in SDs.
        double mz_delta; //!< MZ difference between peaks.
        double mz_lower; //!< MZ lower bound.
        double mz_upper; //!< MZ upper bound.
        double confidence; //!< Confidence for keeping score.  In Standard Deviations.
        int num_threads;
        int input_spectrum_cache_size; //!< Size of input spectrum cache in number of spectra. 
        std::string in_file; //!< Path to input file.
        std::string out_file; //!< Path to output file.

        Options(int argc, char *argv[]);
};

#endif
