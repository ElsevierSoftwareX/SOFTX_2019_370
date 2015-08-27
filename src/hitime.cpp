#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include "options.h"
#include "constants.h"
#include "vector.h"
#include "score.h"

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

   // Calculate number of spectra for each window
   int half_window = ceil(opts.rt_sigma * opts.rt_width / std_dev_in_fwhm);

   // Setup output stream
   std::ofstream outfile;
   outfile.open(opts.out_file);
   outfile.precision(12);

   double_2d score;

   for (Size n = 0; n < map.getNrSpectra(); n++)
   {
       // Score the window
       score = score_spectra(map, n, half_window, opts);

       MSSpectrum<> center_spectrum = map.getSpectrum(n);

       // Output the results
       write_scores(score, center_spectrum, outfile, opts);
   }

   // Close ouput stream
   outfile.close();
   return 0;
}
