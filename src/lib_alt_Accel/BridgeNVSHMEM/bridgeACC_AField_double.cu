/*!
      @file    bridgeACC_AField_double.cu
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2225 $
*/

#include "inline/define_params.h"
#include "inline/define_index.h"
#include "bridgeACC.h"

#include "bridgeACC_AField.h"


typedef double real_t;


namespace BridgeACC {

#include "src/afield_cuda-inc.h"
#include "src/aindex_eo_alt_cuda-inc.h"
#include "src/afield_dd_cuda-inc.h"
#include "src/afield_Gauge_cuda-inc.h"
#include "src/shiftAField_lex_cuda-inc.h"

}

//============================================================END=====
