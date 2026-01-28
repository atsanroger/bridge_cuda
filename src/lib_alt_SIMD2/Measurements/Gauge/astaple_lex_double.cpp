/*!
      @file    astaple_lex_double.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib_alt_SIMD2/Measurements/Gauge/astaple_lex.h"

#include <string>

#include "lib/ResourceManager/threadManager_OpenMP.h"

#include "lib_alt_SIMD2/Field/index_lex_alt.h"
#include "lib_alt_SIMD2/Field/afield.h"

#define  VLEN    VLEND
#define  VLEN2   VLEND2

// inline function files
#include "lib_alt_SIMD2/inline/afield_th-inc.h"

#include "lib_alt_SIMD2/inline/vsimd_double-inc.h"
#include "lib_alt_SIMD2/inline/vsimd_common_double-inc.h"

#if defined USE_GROUP_SU3
#include "lib_alt_SIMD2/inline/define_params_SU3.h"
#include "lib_alt_SIMD2/inline/vsimd_common_SU3_double-inc.h"
#endif

// function template files
#include "lib_alt_SIMD2/Field/afield-inc.h"
#include "lib_alt_SIMD2/Field/afield_Gauge-inc.h"

#include "lib_alt_SIMD2/Measurements/Gauge/astaple_lex-tmpl.h"

template<>
const std::string AStaple_lex<AField<double,SIMD2> >::class_name
                              = "AStaple_lex<AField<double,SIMD2> >";

// explicit instanciation.
template class AStaple_lex<AField<double,SIMD2> >;


//============================================================END=====
