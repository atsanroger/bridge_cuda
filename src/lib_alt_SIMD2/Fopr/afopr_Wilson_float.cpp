/*!
      @file    afopr_Wilson_float.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib_alt_SIMD2/Fopr/afopr_Wilson.h"

#include "lib/ResourceManager/threadManager_OpenMP.h"

#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/index_lex_alt.h"

#include "lib_alt_SIMD2/inline/define_params.h"

#define  VLEN    VLENS
#define  VLEN2   VLENS2

typedef float real_t;

// inline function files
#include "lib_alt_SIMD2/inline/afield_th-inc.h"

#include "lib_alt_SIMD2/inline/vsimd_float-inc.h"
#include "lib_alt_SIMD2/inline/vsimd_common_float-inc.h"

#if defined USE_GROUP_SU3
#include "lib_alt_SIMD2/inline/vsimd_Wilson_SU3_float-inc.h"
#endif

// function template files
#include "lib_alt_SIMD2/Field/afield-inc.h"

// C-code
#include "lib_alt_SIMD2/Fopr/mult_Wilson_simd2-inc.h"

// template class
#include "lib_alt_SIMD2/Fopr/afopr_Wilson-tmpl.h"

template<>
const std::string AFopr_Wilson<AField<float,SIMD2> >::class_name
                            = "AFopr_Wilson<AField<float,SIMD2> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AFopr<AField<float,SIMD2> >::Factory_params::Register(
                               "Wilson", create_object_with_params);
}
#endif


// explicit instanciation.
template class AFopr_Wilson<AField<float,SIMD2> >;

//============================================================END=====
