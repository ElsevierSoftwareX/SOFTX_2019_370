#include "options.h"
#include "score.h"

int main(int argc, char** argv)
{
   Options opts(argc, argv);
   Scorer scorer(opts.debug, opts.intensity_ratio, opts.rt_width,
      opts.mz_width, opts.mz_delta, opts.confidence,
      opts.num_threads, opts.input_spectrum_cache_size, opts.in_file, opts.out_file);
   return 0;
}
