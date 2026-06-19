/*!
      @file    bridgeQXS_Clover_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

#define  VLEN     VLEND
#define  VLENX    VLENXD
#define  VLENY    VLENYD

typedef double real_t;

#include "lib_alt_QXS/inline/vsimd_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_Wilson_SU3_double-inc.h"

#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Clover.h"

#include "src/mult_common_parts_qxs-inc.h"
#include "src/mult_Wilson_parts_qxs-inc.h"
#include "src/mult_Clover_parts_qxs-inc.h"

#include "src/mult_Clover_qxs-inc.h"
#include "src/mult_Clover_dd_qxs-inc.h"
