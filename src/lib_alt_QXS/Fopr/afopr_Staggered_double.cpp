/*!
      @file    afopr_Staggered_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#include "lib_alt_QXS/Fopr/afopr_Staggered.h"

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

#define  VLEN     VLEND
#define  VLENX    VLENXD
#define  VLENY    VLENYD

typedef double real_t;

// inline function files
#include "lib_alt_QXS/inline/vsimd_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_double-inc.h"

#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_Gauge-inc.h"

#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Staggered.h"

#include "lib_alt_QXS/Fopr/mult_common_th-inc.h"


// template definition
#include "lib_alt_QXS/Fopr/afopr_Staggered-tmpl.h"

template<>
const std::string AFopr_Staggered<AField<double, QXS> >::class_name
  = "AFopr_Staggered<AField<double,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AFopr<AField<double, QXS> >::Factory_params::Register(
    "Staggered", create_object_with_params);
}
#endif

// explicit instanciation
template class AFopr_Staggered<AField<double, QXS> >;

//============================================================END=====
