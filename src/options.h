#include <unistd.h>

//! @brief Print program usage information.
void show_usage(char *cmd);

class Options {

    public:
        float intensity_ratio; //!< Intensity ratio between lo and hi peaks.
        float rt_width; //!< Retention time FWHM in scans.
        float rt_sigma; //!< Boundary for RT width in SDs.
        float ppm; //!< MZ tolerance in PPM.
        float mz_width; //!< MZ FWHM in PPM.
        float mz_sigma; //!< Boundary for MZ in SDs.
        float mz_delta; //!< MZ difference between peaks.
        float min_sample; //!< Minimum number of points required in each region.
        bool full_out; //!< Output all points (including zero scores).
        std::string mzML_file; //!< Path to input file.
        std::string out_file; //!< Path to output file.

        Options(int argc, char *argv[]);
};


