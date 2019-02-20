#ifndef TAGE_BASE_PREDICTOR
#define TAGE_BASE_PREDICTOR

#include "simple_bimodal_table.h"

class TageBasePredictor
   : public SimpleBimodalTable
{

public:
   TageBasePredictor()
      : SimpleBimodalTable(4096)
   {}

};

#endif /* TAGE_BASE_PREDICTOR */
