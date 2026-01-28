
/*!
      @file    bridgeACC_float.cu
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2225 $
*/

#include "inline/define_params.h"
#include "inline/define_index.h"

#include "bridgeACC.h"

typedef float real_t;

namespace BridgeACC {

#include "src/bridge_cuda-inc.h"

}

//============================================================END=====
