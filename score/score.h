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
  unsigned int num_spectra;
  unsigned int num_threads;
  int current_spectrum;
  int half_window;
  Options *opts;
  OnDiscPeakMap input_map;
  MSExperiment output_map;
  InputSpectrumCache input_spectrum_cache;
  int get_next_spectrum_todo(void);
  void write_spectrum(MSSpectrum<> spectrum);
  MSSpectrum<> read_spectrum(int spectrum_id);
  double_vect score_spectra(int centre_idx);

public:
  Scorer(Options *o);
  void score_worker(int thread_count);
};
#endif
