/*!
      @file    bridgeACC_Domainwall_5din_float.cu
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "inline/define_params.h"
#include "inline/define_index.h"

#include "bridgeACC.h"

#include "bridgeACC_AField.h"


typedef float real_t;
#include "inline/mult_Wilson_cuda_inline-inc.h"
//#include "src/mult_Wilson_cuda2-inc.h"

namespace BridgeACC {

#define SU3_3RD_ROW_RECONST

#include "src/mult_Domainwall_5din_cuda-inc.h"
#include "src/mult_Domainwall_5din_eo_cuda-inc.h"
#include "src/mult_Domainwall_5din_dd_cuda-inc.h"
#include "src/mult_Domainwall_5din_cuda2-inc.h"

}

//============================================================END=====
