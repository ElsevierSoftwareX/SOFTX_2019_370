#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <thread>
#include "options.h"
#include "constants.h"
#include "vector.h"
#include "score.h"

#include <iostream>
#include <iterator>
using namespace std;
using namespace OpenMS;

int main(int argc, char** argv)
{
   // Read user options
   Options opts(argc, argv);

   cout << "I'm here: " << endl;

   IndexedMzMLFileLoader imzml;

   // load data from an indexed MzML file
   OnDiscPeakMap map;
   imzml.load(opts.in_file, map);

   PlainMSDataWritingConsumer* consumer = new PlainMSDataWritingConsumer(opts.out_file);
   consumer->setExpectedSize(map.getNrSpectra(), map.getNrChromatograms());
   //consumer->setExperimentalSettings(map.getExperimentalSettings());

   
   MSSpectrum<> s1 = map.getSpectrum(0);
   MSSpectrum<> s2 = map.getSpectrum(1);
   MSSpectrum<> s3 = map.getSpectrum(2);

   consumer->consumeSpectrum(s1);
   consumer->consumeSpectrum(s2);
   consumer->consumeSpectrum(s3);
   delete consumer;


   //IndexedMzMLFileLoader imzml;
   //MzMLFile mzml;

   // load data from an indexed MzML file
   //MSExperiment input_map;
   //MSExperiment output_map;

   //imzml.load(opts.in_file, input_map);
   //mzml.load(opts.in_file, input_map);


   return 0;
}
