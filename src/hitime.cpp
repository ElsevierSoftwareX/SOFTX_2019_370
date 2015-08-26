#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include "options.h"

using namespace OpenMS;
using namespace std;

int main(int argc, char** argv)
{
   // Read user options
   Options opts(argc, argv);

   IndexedMzMLFileLoader imzml;
   //MzMLFile mzml;

   // load data from an indexed MzML file
   OnDiscMSExperiment<> map;
   //MSExperiment<> map;

   imzml.load(opts.mzML_file, map);
   //mzml.load(opts.mzML_file, map);

   // Get the first spectrum in memory, do some constant (non-changing) data processing
   MSSpectrum<> s = map.getSpectrum(0);
   std::cout << "There are " << map.getNrSpectra() << " spectra in the input file." << std::endl;
   std::cout << "The first spectrum has " << s.size() << " peaks." << std::endl;

   // store the (unmodified) data in a different file
   //imzml.store("Tutorial_FileIO_output.mzML", map);

   // Calculate number of spectra for each window
   // XXX fix magic number
   int half_window = ceil(opts.rt_sigma * opts.rt_width / 2.355);


/*
   // Access the input file using pwiz
   pwiz::msdata::FullReaderList readers;
   pwiz::msdata::MSDataFile msd(opts.mzML_file, &readers);
   pwiz::msdata::SpectrumList& spectrumList = *msd.run.spectrumListPtr;

   int rt_len = spectrumList.size();
   double_2d score;
   pwiz::msdata::SpectrumPtr centre_vect;

   // Setup output stream
   std::ofstream outfile;
   outfile.open(opts.out_file);
   outfile.precision(12);

   // Loop over the spectra
   for (int centre_rt = 0; centre_rt < rt_len; ++centre_rt) {

       // Select the centre spectrum
       centre_vect = spectrumList.spectrum(centre_rt, opts.getBinaryData);

       // Score the window
       score = score_spectra(msd, centre_rt, half_window, opts);

       // Output the results
       write_scores(score, centre_vect, outfile, opts);
   }

   // Close ouput stream
   outfile.close();
*/

   std::cout << "Done!" << std::endl;
   return 0;

}

/*! Scores are output in CSV format with the following fields: retention time,
 * mz, intensity, minimum score, correlAB, correlA0, correlB0, correl1r.
 *
 * @param scores 2D vector of scores returned by score_spectra.
 * @param raw_data Spectrum pointer to the raw central vector.
 * @param out_stream Stream to write output too.
 * @param opts User defined Options object.
 */
/*
void write_scores(double_2d scores, pwiz::msdata::SpectrumPtr raw_data,
                  std::ofstream& out_stream, Options opts)
{
    // Get central spectrum retention time
    pwiz::msdata::SpectrumInfo spectrum_info;
    spectrum_info.update(*raw_data, opts.getBinaryData);
    double rt = spectrum_info.retentionTime;

    // Get raw MZ/Intensity pairs
    std::vector<pwiz::msdata::MZIntensityPair> raw_pairs;
    raw_data->getMZIntensityPairs(raw_pairs);

    // Write output
    for (size_t idx = 0; idx < raw_pairs.size(); ++idx) {
        double mz  = raw_pairs[idx].mz;
        double amp = raw_pairs[idx].intensity;
        double ms  = scores[0][idx];
        double AB  = scores[1][idx];
        double A0  = scores[2][idx];
        double B0  = scores[3][idx];
        double r1  = scores[4][idx];

        if (opts.full_out == true) {
            out_stream << rt << ", " << mz << ", " << amp << ", "
                       << ms << ", " << AB << ", " << A0 << ", "
                       << B0 << ", " << r1 << std::endl;
        } else {
            if (ms > 0.0) {
                out_stream << rt << ", " << mz << ", " << amp << ", "
                           << ms << ", " << AB << ", " << A0 << ", "
                           << B0 << ", " << r1 << std::endl;
            }
        }
    }
}
*/
