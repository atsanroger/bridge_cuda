/*!
        @file    aforce_F_CloverTerm_double.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#include "lib_alt_SIMD2/Force/Fermion/aforce_F_CloverTerm.h"

#include<cassert>

#include "lib/ResourceManager/threadManager_OpenMP.h"
#include "lib/IO/bridgeIO.h"
#include "lib/Tools/gammaMatrixSet.h"
using Bridge::vout;

// include files in alt-code dorectories
#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/index_lex_alt.h"
#include "lib_alt_SIMD2/Field/shiftAField_lex.h"
#include "lib_alt_SIMD2/Measurements/Gauge/astaple_lex.h"


#define  VLEN    VLEND
#define  VLEN2   VLEND2

#include "lib_alt_SIMD2/inline/afield_th-inc.h"

#include "lib_alt_SIMD2/inline/vsimd_double-inc.h"
#include "lib_alt_SIMD2/inline/vsimd_common_double-inc.h"

#if defined USE_GROUP_SU3
#include "lib_alt_SIMD2/inline/vsimd_Wilson_SU3_double-inc.h"
#include "lib_alt_SIMD2/inline/vsimd_common_SU3_double-inc.h"
#endif

#include "lib_alt_SIMD2/Fopr/mult_common_th_simd2-inc.h"
#include "lib_alt_SIMD2/Fopr/mult_Wilson_simd2_parts-inc.h"

#include "lib_alt_SIMD2/Field/afield_4spinor-inc.h"
#include "lib_alt_SIMD2/Field/afield_Gauge-inc.h"

// function template files
#include "lib_alt_SIMD2/Field/afield-inc.h"

#include "lib_alt_SIMD2/Force/Fermion/aforce_F_CloverTerm-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,SIMD2>.
template<>
const std::string AForce_F_CloverTerm<AField<double,SIMD2> >::class_name
= "AForce_F_CloverTerm<Afield<double,SIMD2> >";

#ifdef USE_FACTORY_AUTOREGISTER
// This class is assumed not to be constructed via Factory.
/*
namespace {
  bool init1 = AForce_F_CloverTerm<AField_dev<double,SIMD2> >::register_factory();
}
*/
#endif

template class AForce_F_CloverTerm<AField<double,SIMD2> >;

//============================================================END=====
