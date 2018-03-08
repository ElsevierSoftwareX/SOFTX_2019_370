#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include <thread>
#include "options.h"
#include "constants.h"
#include "vector.h"
#include "maxima.h"

#include <iostream>
#include <iterator>
using namespace std;
using namespace OpenMS;

int main(int argc, char** argv)
{
   // Read user options
   Options opts(argc, argv);

   //IndexedMzMLFileLoader imzml;
   MzMLFile mzml;

   // load data from an indexed MzML file
   MSExperiment input_map;
   MSExperiment output_map;

   //imzml.load(opts.in_file, input_map);
   mzml.load(opts.in_file, input_map);

   // Calculate number of spectra for each window
   int half_window = ceil(opts.rt_sigma * opts.rt_width / std_dev_in_fwhm);

   int num_spectra = input_map.getNrSpectra();
   int spectra_per_thread = num_spectra / opts.num_threads;
   vector<thread> threads(opts.num_threads);

   int low_spectrum = 0;
   int high_spectrum = spectra_per_thread;

   cout.precision(17);

   //cout << "Num threads: " << opts.num_threads << endl;
   //cout << "Spectra per thread: " << spectra_per_thread << endl;

   for (int thread_count = 0; thread_count < opts.num_threads; thread_count++)
   {
      //cout << thread_count << " " << low_spectrum << " " << high_spectrum << endl;
      threads[thread_count] = thread {score_worker, std::ref(input_map), std::ref(output_map), half_window, opts, low_spectrum, high_spectrum};
      // score_worker (input_map, output_map, half_window, opts, low_spectrum, high_spectrum);
      low_spectrum += spectra_per_thread;

      // The high spectrum is special for the last thread.
      if (thread_count == opts.num_threads - 1)
      {
          high_spectrum = num_spectra - 1;
      }
      else
      {
          high_spectrum += spectra_per_thread;
      }
   }

   // Wait for all the threads to complete
   for (int thread_count = 0; thread_count < opts.num_threads; thread_count++)
   {
       threads[thread_count].join();
   }

   // Write the output scored mzml file
   //imzml.store(opts.out_file, output_map);
   mzml.store(opts.out_file, output_map);

   // Close ouput stream
   // outfile.close();

   return 0;
}
