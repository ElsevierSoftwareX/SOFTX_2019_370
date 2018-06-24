#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include <thread>
#include "options.h"
#include "constants.h"
#include "vector.h"
#include "score.h"

#include <map>
#include <iostream>
#include <iterator>
using namespace std;
using namespace OpenMS;

int main(int argc, char** argv)
{
   // Read user options
   Options opts(argc, argv);

   IndexedMzMLFileLoader mzml;

   // load data from an indexed MzML file
   OnDiscPeakMap input_map;
   MSExperiment output_map;

   InputSpectrumCache input_spctrum_cache;

   mzml.load(opts.in_file, input_map);

   int half_window = ceil(opts.rt_sigma * opts.rt_width / std_dev_in_fwhm);

   int num_spectra = input_map.getNrSpectra();
   vector<thread> threads(opts.num_threads);

   int current_spectrum = 0;

   cout << "Num threads: " << opts.num_threads << endl;
   cout << "Num spectra: " << num_spectra << endl;

   for (int thread_count = 0; thread_count < opts.num_threads; thread_count++)
   {
      threads[thread_count] = thread {score_worker, std::ref(input_map), std::ref(output_map), std::ref(input_spctrum_cache), half_window, opts, &current_spectrum, num_spectra, thread_count};
   }

   // Wait for all the threads to complete
   for (int thread_count = 0; thread_count < opts.num_threads; thread_count++)
   {
       threads[thread_count].join();
   }

   // Write the output scored mzml file
   mzml.store(opts.out_file, output_map);

   return 0;
}
