#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include "options.h"
#include "constants.h"
#include "vector.h"
#include "score.h"

#include <iostream>
#include <iterator>
using namespace std;
using namespace OpenMS;

int main(int argc, char** argv)
{
   // Read user options
   Options opts(argc, argv);

   IndexedMzMLFileLoader imzml;

   // load data from an indexed MzML file
   OnDiscMSExperiment<> input_map;
   MSExperiment<> output_map;

   imzml.load(opts.in_file, input_map);

   // Calculate number of spectra for each window
   int half_window = ceil(opts.rt_sigma * opts.rt_width / std_dev_in_fwhm);

   double_2d score;

   for (Size n = 0; n < input_map.getNrSpectra(); n++)
   {
       // Score the window
       score = score_spectra(input_map, n, half_window, opts);

       MSSpectrum<> input_spectrum = input_map.getSpectrum(n);
       // Copy the input spectrum to the output spectrum
       // XXX this must be a deep copy
       MSSpectrum<> output_spectrum = MSSpectrum<Peak1D>(input_spectrum);

       // Write the score for each peak into the output spectrum (mz is kept the same in the output).
       for (int index = 0; index < input_spectrum.size(); ++index)
       {
           output_spectrum[index].setIntensity(score[0][index]);
       }
       output_map.addSpectrum(output_spectrum);

       // Output the results
       // write_scores(score, input_spectrum, outfile);
   }

   // Write the output scored mzml file
   imzml.store(opts.out_file, output_map);

   // Close ouput stream
   // outfile.close();

   return 0;
}
