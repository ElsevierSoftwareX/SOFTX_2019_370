#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
using namespace OpenMS;

//! @brief Calculate correlation scores for a window of data.
double_2d score_spectra(OnDiscMSExperiment<> map, int centre_idx,
                        int half_window, Options opts);

//! @brief Write scores to output file
void write_scores(double_2d scores, MSSpectrum<> raw_data, std::ofstream& out_stream);
