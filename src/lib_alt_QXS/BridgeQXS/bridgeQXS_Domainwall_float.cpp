/*!
      @file    bridgeQXS_Domainwall_float.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2581 $
*/

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

typedef float real_t;

#include "lib_alt_QXS/inline/vsimd_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_Wilson_SU3_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_Domainwall_SU3_float-inc.h"

#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Domainwall.h"

#include "src/mult_common_parts_qxs-inc.h"
#include "src/mult_Wilson_parts_qxs-inc.h"
#include "src/mult_Wilson_parts_qxs2-inc.h"
#include "src/mult_Wilson_eo_parts_qxs-inc.h"

#include "src/mult_Domainwall_5din_eo_qxs-inc.h"
#include "src/mult_Domainwall_5din_qxs-inc.h"
#include "src/mult_Domainwall_5din_dd_qxs-inc.h"
