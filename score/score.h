#ifndef HITIME_SCORE_H
#define HITIME_SCORE_H

#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include "options.h"
#include "vector.h"
using namespace OpenMS;
using namespace std;

typedef std::map<int, MSSpectrum<>> InputSpectrumCache;

class Scorer 
{
private:
   // attributes
   unsigned int num_spectra;
   int current_spectrum;
   int half_window;
   bool debug;
   double intensity_ratio;
   double rt_width;
   double rt_sigma;
   double ppm;
   double mz_width;
   double mz_sigma;
   double mz_delta;
   double min_sample;
   unsigned int num_threads;
   string in_file;
   string out_file;
   OnDiscPeakMap input_map;
   MSExperiment output_map;
   InputSpectrumCache input_spectrum_cache;

   // methods
   int get_next_spectrum_todo(void);
   void write_spectrum(MSSpectrum<> spectrum);
   MSSpectrum<> read_spectrum(int spectrum_id);
   double_vect score_spectra(int centre_idx);

public:
   Scorer(bool debug, double intensity_ratio, double rt_width, double rt_sigma,
         double ppm, double mz_width, double mz_sigma, double mz_delta,
         double min_sample, int num_threads, string in_file, string out_file);
  void score_worker(int thread_count);
};
#endif
