/*!
      @file    bridgeACC_Staggered_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "define_params.h"
#include "define_index.h"

typedef double real_t;

namespace BridgeACC {

#include "src/mult_Staggered_openacc2R-inc.h"
#include "src/mult_Staggered_eo_openacc2R-inc.h"
#include "src/mult_Staggered_openacc-inc.h"
#include "src/mult_Staggered_eo_openacc-inc.h"

}

//============================================================END=====
