#include "options.h"
#include "score.h"

int main(int argc, char** argv)
{
   Options opts(argc, argv);
   Scorer scorer(opts.debug, opts.intensity_ratio, opts.rt_width, opts.rt_sigma,
      opts.ppm, opts.mz_width, opts.mz_sigma, opts.mz_delta, opts.min_sample,
      opts.num_threads, opts.input_spectrum_cache_size, opts.in_file, opts.out_file);
   return 0;
}
