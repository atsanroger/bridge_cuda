/*!
      @file    afopr_Clover_eo_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib_alt_SIMD2/Fopr/afopr_Clover_eo.h"

#include "lib/ResourceManager/threadManager_OpenMP.h"

#include "lib_alt_SIMD2/inline/define_params.h"

#define  VLEN   VLEND
#define  VLEN2  VLEND2

#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/index_lex_alt.h"
#include "lib_alt_SIMD2/Field/index_eo_alt.h"

// inline function files
#include "lib_alt_SIMD2/inline/afield_th-inc.h"

#include "lib_alt_SIMD2/inline/vsimd_double-inc.h"
#include "lib_alt_SIMD2/inline/vsimd_common_double-inc.h"

#if defined USE_GROUP_SU3
#include "lib_alt_SIMD2/inline/vsimd_Wilson_SU3_double-inc.h"
#endif

#include "lib_alt_SIMD2/Fopr/mult_common_th_simd2-inc.h"
#include "lib_alt_SIMD2/Fopr/mult_Wilson_simd2_parts-inc.h"

// function template files
#include "lib_alt_SIMD2/Field/afield-inc.h"
#include "lib_alt_SIMD2/Field/index_eo_alt-inc.h"

#include "lib_alt_SIMD2/Fopr/afopr_Clover_eo-tmpl.h"


template<>
const std::string AFopr_Clover_eo<AField<double,SIMD2> >::class_name
                           = "AFopr_Clover_eo<AField<double,SIMD2> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AFopr<AField<double,SIMD2> >::Factory_params::Register(
                          "Clover_eo", create_object_with_params);
}
#endif

// explicit instanciation.
template class AFopr_Clover_eo<AField<double,SIMD2> >;

//============================================================END=====
