/*!
      @file    bridgeQXS_Staggered_double.cpp
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

// inline functions
#include "lib_alt_QXS/inline/vsimd_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_Wilson_SU3_double-inc.h"

#include "src/mult_common_parts_qxs-inc.h"
#include "src/mult_Staggered_parts_qxs-inc.h"

// forward declaration
#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Staggered.h"

// instantiation
#include "src/mult_Staggered_qxs-inc.h"
#include "src/mult_Staggered_eo_qxs-inc.h"
