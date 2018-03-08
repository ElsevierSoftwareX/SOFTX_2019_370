#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
using namespace OpenMS;

void score_worker(MSExperiment &input_map, MSExperiment &output_map, int half_window, Options opts, int low_spectra, int high_spectra);

double_vect score_spectra(MSExperiment &map, int centre_idx, int half_window, Options opts);
