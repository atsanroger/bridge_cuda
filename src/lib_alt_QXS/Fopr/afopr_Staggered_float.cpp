/*!
      @file    afopr_Staggered_float.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#include "lib_alt_QXS/Fopr/afopr_Staggered.h"

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

typedef float real_t;

// inline function files
#include "lib_alt_QXS/inline/vsimd_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_float-inc.h"

#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_Gauge-inc.h"

#include "lib_alt_QXS/Fopr/mult_common_th-inc.h"

#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Staggered.h"


// template definition
#include "lib_alt_QXS/Fopr/afopr_Staggered-tmpl.h"

template<>
const std::string AFopr_Staggered<AField<float, QXS> >::class_name
  = "AFopr_Staggered<AField<float,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AFopr<AField<float, QXS> >::Factory_params::Register(
    "Staggered", create_object_with_params);
}
#endif

// explicit instanciation
template class AFopr_Staggered<AField<float, QXS> >;

//============================================================END=====
