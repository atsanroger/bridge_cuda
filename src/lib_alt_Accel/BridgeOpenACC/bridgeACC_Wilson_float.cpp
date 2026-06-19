/*!
      @file    bridgeACC_Wilson_float.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2574 $
*/

#include <vector>

#include "bridge_defs.h"

#include "define_params.h"
#include "define_index.h"

typedef float real_t;

namespace BridgeACC {

#define SU3_3RD_ROW_RECONST

#include "bridgeACC_Wilson.h"

#include "src/mult_Wilson_openacc2-inc.h"
#include "src/mult_Wilson_eo_openacc2-inc.h"
#include "src/mult_Wilson_openacc-inc.h"
#include "src/mult_Wilson_eo_openacc-inc.h"
#include "src/mult_Wilson_dd_openacc-inc.h"

}

//============================================================END=====
