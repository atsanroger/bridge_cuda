/*!
      @file    afopr_Wilson_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#include "lib_alt_QXS/Fopr/afopr_Wilson.h"

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

#define  VLEN     VLEND
#define  VLENX    VLENXD
#define  VLENY    VLENYD

typedef double real_t;

#include "lib_alt_QXS/inline/vsimd_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_Wilson_SU3_double-inc.h"

#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_Gauge-inc.h"

#include "lib_alt_QXS/Fopr/mult_common_th-inc.h"
#include "lib_alt_QXS/Fopr/mult_Wilson_parts_qxs_org-inc.h"

#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Wilson.h"


// template definition
#include "lib_alt_QXS/Fopr/afopr_Wilson-tmpl.h"

template<>
const std::string AFopr_Wilson<AField<double, QXS> >::class_name
  = "AFopr_Wilson<AField<double,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AFopr<AField<double, QXS> >::Factory_params::Register(
    "Wilson", create_object_with_params);
}
#endif

// explicit instanciation
template class AFopr_Wilson<AField<double, QXS> >;

//============================================================END=====
