/*!
      @file    bridgeACC_Staggered_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2225 $
*/

#include "define_params.h"
#include "define_index.h"

typedef float real_t;

namespace BridgeACC {

#include "src/afield_openacc-inc.h"
#include "src/afield_Gauge_openacc-inc.h"
#include "src/index_eo_alt_openacc-inc.h"
#include "src/shiftAField_lex_openacc-inc.h"
#include "src/afield_dd_openacc-inc.h"

  // copy with double/float conversion
#include "src/afield_copy_openacc-inc.h"

}

//============================================================END=====
