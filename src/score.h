/*
#include <iostream>
#include <unistd.h>
#include <math.h>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <functional>
#include <fstream>
#include <numeric>
#include <iostream>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>

using namespace std;
*/

#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include "vector.h"
using namespace OpenMS;

//! @brief Calculate correlation scores for a window of data.

/*
double_2d score_spectra(pwiz::msdata::MSDataFile &msd, int centre_idx,
                        int half_window, Options opts);
*/

double_2d score_spectra(OnDiscMSExperiment<> map, int centre_idx,
                        int half_window, Options opts);
