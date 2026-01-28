/*!
      @file    bridgeACC_Domainwall_coarse_float.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2578 $
*/

#include "define_params.h"
#include "define_index.h"

typedef float real_t;

namespace BridgeACC {

#include "src/mult_Domainwall_coarse_openacc-inc.h"

}

//============================================================END=====
