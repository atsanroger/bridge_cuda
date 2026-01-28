/*!
      @file    bridgeACC_Domainwall_coarse_double.cu
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2578 $
*/

#include "inline/define_params.h"
#include "inline/define_index.h"

#include "bridgeACC.h"
#include "bridgeACC_AField.h"

using namespace std;

typedef double real_t;

#include "inline/mult_Wilson_cuda_inline-inc.h"

namespace BridgeACC {

#include "src/mult_Domainwall_coarse_cuda-inc.h"

}

//============================================================END=====
