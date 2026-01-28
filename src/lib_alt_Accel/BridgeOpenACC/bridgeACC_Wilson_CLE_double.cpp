/*!
      @file    bridgeACC_Wilson_CLE_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "define_params.h"
#include "define_index.h"

typedef double real_t;

namespace BridgeACC {

//#define SU3_3RD_ROW_RECONST
// For CLE, do not reconstruct !

#include "src/mult_Wilson_CLE_openacc2-inc.h"

}

//============================================================END=====
