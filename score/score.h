#ifndef HITIME_SCORE_H
#define HITIME_SCORE_H

#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <queue>
#include "options.h"
#include "vector.h"
#include "lru_cache.h"

using namespace OpenMS;
using namespace std;

struct IndexSpectrumOrder 
{
   bool operator()(pair<int, PeakSpectrum> const& a, pair<int, PeakSpectrum> const& b) const
   {
       return a.first > b.first;
   }
};

typedef shared_ptr<PeakSpectrum> PeakSpectrumPtr;
typedef cache::lru_cache<Int, PeakSpectrumPtr> SpectrumLRUCache;
typedef pair<int, PeakSpectrum> IndexSpectrum;
typedef priority_queue<IndexSpectrum, vector<IndexSpectrum>, IndexSpectrumOrder> SpectrumQueue;

class Scorer 
{
private:
   // attributes
   unsigned int num_spectra;
   int current_spectrum_id;
   int next_output_spectrum_id;
   int half_window;
   OpenMS::Size local_rows;
   bool debug;
   double intensity_ratio;
   double rt_width;
   double rt_sigma;
   double mz_width;
   double mz_sigma;
   double mz_delta;
   double mz_upper;
   double mz_lower;
   double min_sample;
   double confidence;
   unsigned int num_threads;
   string in_file;
   string out_file;
   OnDiscPeakMap input_map;
   PlainMSDataWritingConsumer spectrum_writer;
   SpectrumLRUCache input_spectrum_cache;
   SpectrumQueue output_spectrum_queue;
   std::ofstream csv_fs;
   
   // methods
   int get_next_spectrum_todo(void);
   void put_spectrum(int spectrum_id, PeakSpectrum spectrum);
   PeakSpectrumPtr get_spectrum(int spectrum_id);
   PeakSpectrum score_spectra(int centre_idx);
   PeakSpectrum local_max_spectra(int centre_idx);
   void collect_local_rows(int, double_2d&, double_2d&);
   void collect_window_data(double,
                  double_vect&, double, double, double_2d&, double_2d&,
                  double, double, double_vect&, double_vect&);
   bool local_max_data(double,
                  double_2d&, double_2d&,
                  double, double);

public:
   Scorer(bool debug, double intensity_ratio, double rt_width, 
         double mz_width, double mz_delta, double mz_lower, double mz_upper, double confidence,
         int num_threads, int input_spectrum_cache_size,
         string in_file, string out_file);
  void score_worker(int thread_count);
};
#endif
