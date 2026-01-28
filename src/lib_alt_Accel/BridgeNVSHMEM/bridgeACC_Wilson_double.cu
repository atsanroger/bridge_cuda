/*!
      @file    bridgeACC_Wilson_double.cu
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "inline/define_params.h"
#include "inline/define_index.h"
#include "bridgeACC.h"

#include "bridgeACC_AField.h"


typedef double real_t;


namespace BridgeACC {

#define SU3_3RD_ROW_RECONST

#include "src/mult_Wilson_cuda-inc.h"
#include "src/mult_Wilson_cuda2-inc.h"
#include "src/mult_Wilson_dd_cuda-inc.h"


}

//============================================================END=====
