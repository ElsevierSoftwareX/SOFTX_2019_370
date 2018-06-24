#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
using namespace OpenMS;

typedef std::map<int, MSSpectrum<>> InputSpectrumCache;

void score_worker(OnDiscPeakMap &input_map, MSExperiment &output_map, InputSpectrumCache &input_spectrum_cache, int half_window, Options opts, int *current_spectra, int num_spectra, int thread_count);

double_vect score_spectra(OnDiscPeakMap &map, int centre_idx, int half_window, Options opts);
