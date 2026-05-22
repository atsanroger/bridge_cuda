/*!
      @file    bridgeACC_Domainwall_5din_double.cu
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
#include "inline/mult_Wilson_cuda_inline-inc.h"
//#include "src/mult_Wilson_cuda2-inc.h"

namespace BridgeACC {

#define SU3_3RD_ROW_RECONST
//#define USE_5D_HOPB_V2  // shared-memory gauge link sharing for 5D EO hopb kernel
//#define USE_LUINV_V2    // thread-local vt_fwd[NS][ND] for LU(dag)inv (drops vp round-trip)

#include "src/mult_Domainwall_5din_cuda-inc.h"
#include "src/mult_Domainwall_5din_eo_cuda-inc.h"
#include "src/mult_Domainwall_5din_dd_cuda-inc.h"
#include "src/mult_Domainwall_5din_cuda2-inc.h"
#include "src/mult_Domainwall_5din_eo_inv_cuda-inc.h"
#include "src/mult_Domainwall_5din_cuda_qdw-inc.h"
#include "src/mult_Domainwall_5din_eo_cuda_qdw-inc.h"
#include "src/mult_Domainwall_5din_eo_inv_cuda_qdw-inc.h"

}

//============================================================END=====
