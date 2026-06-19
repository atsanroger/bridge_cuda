/*!
      @file    bridgeQXS_Clover_coarse_double.cpp
      @brief
      @author  Issaku Kanamori (kanamori)
      @date    $LastChangedDate: 2023-03-20 10:52:44 +0900 (Mon, 20 Mar 2023) $
      @version $LastChangedRevision: 2499 $
*/

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

#define  VLEN     VLEND
#define  VLENX    VLENXD
#define  VLENY    VLENYD

typedef double real_t;

//#include "lib_alt_QXS/inline/define_index.h"
#include "lib_alt_QXS/inline/vsimd_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_Wilson_SU3_double-inc.h"

#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Clover_coarse.h"

#include "src/mult_Clover_coarse_parts_qxs-inc.h"

#include "src/mult_Clover_coarse_qxs-inc.h"
