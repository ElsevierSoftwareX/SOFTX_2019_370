#include <unistd.h>
#include <boost/program_options.hpp>

class Options {

    public:
        double rt_width; //!< Retention time FWHM in scans. keep
        double rt_sigma; //!< Boundary for RT width in SDs. keep
        double mz_width; //!< MZ FWHM in PPM. keep
        double mz_sigma; //!< Boundary for MZ in SDs. keep 
        std::string in_file; //!< Path to input file.
        std::string out_file; //!< Path to output file.
        bool debug;
        int num_threads;

        Options(int argc, char *argv[]);
};
