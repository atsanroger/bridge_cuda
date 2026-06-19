/*!
      @file    bridgeACC_Domainwall_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2575 $
*/

#include "lib_alt_Accel/inline/define_params.h"
#include "lib_alt_Accel/inline/define_index.h"

typedef double real_t;

namespace BridgeACC {

#define SU3_3RD_ROW_RECONST

#include "bridgeACC_Domainwall.h"

#include "src/inc/mult_Wilson_openacc_inline-inc.h"

#include "src/mult_Domainwall_5din_openacc-inc.h"
#include "src/mult_Domainwall_5din_eo_openacc-inc.h"
#include "src/mult_Domainwall_5din_dd_openacc-inc.h"
#include "src/mult_Domainwall_5din_openacc2-inc.h"

#include "src/mult_Domainwall_5din_eo_inv_openacc-inc.h"

}

//============================================================END=====
