#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
using namespace OpenMS;

void score_worker(MSExperiment &input_map, MSExperiment &output_map, int half_window, Options opts, int low_spectra, int high_spectra);

//! @brief Calculate correlation scores for a window of data.
double_vect score_spectra(MSExperiment &map, int centre_idx,
                        int half_window, Options opts);

//! @brief Write scores to output file
void write_scores(double_2d scores, MSSpectrum<> raw_data, std::ofstream& out_stream);
