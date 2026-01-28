/*!
      @file    afopr_Wilson_eo_half.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/
#ifdef USE_QXS_FP16

#include "lib_alt_QXS/Fopr/afopr_Wilson_eo.h"

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

#define  VLEN     VLENH
#define  VLENX    VLENXH
#define  VLENY    VLENYH

typedef half real_t;

#include "lib_alt_QXS/inline/vsimd_half-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_half-inc.h"
#include "lib_alt_QXS/inline/vsimd_Wilson_SU3_half-inc.h"

#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_Gauge-inc.h"

#include "lib_alt_QXS/Fopr/mult_common_th-inc.h"
#include "lib_alt_QXS/Fopr/mult_Wilson_qxs_parts-inc.h"

#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Wilson.h"

// template definition
#include "lib_alt_QXS/Fopr/afopr_Wilson_eo-tmpl.h"

template<>
const std::string AFopr_Wilson_eo<AField<half, QXS> >::class_name
  = "AFopr_Wilson_eo<AField<half,QXS> >";


// explicit instanciation.
template class AFopr_Wilson_eo<AField<half, QXS> >;

#endif
//============================================================END=====
