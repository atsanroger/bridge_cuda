/*!
      @file    afopr_Clover_float.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib_alt_SIMD2/Fopr/afopr_Clover.h"

#include "lib/ResourceManager/threadManager_OpenMP.h"
#include "lib/Measurements/Gauge/staple_lex.h"
#include "lib/Tools/gammaMatrixSet.h"
#include "lib/Field/field_G.h"

#include "lib_alt_SIMD2/inline/define_params.h"

#define  VLEN    VLENS
#define  VLEN2   VLENS2

#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/index_lex_alt.h"
#include "lib_alt_SIMD2/Field/shiftAField_lex.h"
#include "lib_alt_SIMD2/Measurements/Gauge/astaple_lex.h"


// inline function files
#include "lib_alt_SIMD2/inline/afield_th-inc.h"

#include "lib_alt_SIMD2/inline/vsimd_float-inc.h"
#include "lib_alt_SIMD2/inline/vsimd_common_float-inc.h"

#if defined USE_GROUP_SU3
#include "lib_alt_SIMD2/inline/vsimd_common_SU3_float-inc.h"
#include "lib_alt_SIMD2/inline/vsimd_Wilson_SU3_float-inc.h"
#endif

#include "lib_alt_SIMD2/Fopr/mult_common_th_simd2-inc.h"
#include "lib_alt_SIMD2/Fopr/mult_Wilson_simd2_parts-inc.h"

// function template files
#include "lib_alt_SIMD2/Field/afield-inc.h"
#include "lib_alt_SIMD2/Field/afield_Gauge-inc.h"

#include "lib_alt_SIMD2/Fopr/afopr_Clover-tmpl.h"


template<>
const std::string AFopr_Clover<AField<float,SIMD2> >::class_name
                            = "AFopr_Clover<AField<float,SIMD2> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AFopr<AField<float,SIMD2> >::Factory_params::Register(
                              "Clover", create_object_with_params);
}
#endif


// explicit instanciation.
template class AFopr_Clover<AField<float,SIMD2> >;

//============================================================END=====
