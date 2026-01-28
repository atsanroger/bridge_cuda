/*!
      @file    bridgeQXS_Clover_coarse_float.cpp
      @brief
      @author  Issaku Kanamori (kanamori)
      @date    $LastChangedDate: 2023-03-20 10:52:44 +0900 (Mon, 20 Mar 2023) $
      @version $LastChangedRevision: 2499 $
*/

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

typedef float real_t;

//#include "lib_alt_QXS/inline/define_index.h"
#include "lib_alt_QXS/inline/vsimd_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_Wilson_SU3_float-inc.h"

#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Clover_coarse.h"

#include "src/mult_Clover_coarse_parts_qxs-inc.h"

#include "src/mult_Clover_coarse_qxs-inc.h"
