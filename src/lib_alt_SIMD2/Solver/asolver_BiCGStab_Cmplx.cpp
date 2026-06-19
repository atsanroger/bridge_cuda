/*!
      @file    asolver_BiCGStab_Cmplx.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib_alt/Solver/asolver_BiCGStab_Cmplx.h"

#include "lib/ResourceManager/threadManager_OpenMP.h"

#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/afield-inc.h"

#include "lib_alt/Solver/asolver_BiCGStab_Cmplx-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,SIMD2>.
template<>
const std::string ASolver_BiCGStab_Cmplx<AField<double,SIMD2> >::class_name
                     = "ASolver_BiCGStab_Cmplx<Afield<double,SIMD2> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_BiCGStab_Cmplx<AField_dev<double,SIMD2> >::register_factory();
}
#endif

template class ASolver_BiCGStab_Cmplx<AField<double,SIMD2> >;

//====================================================================
// explicit instanciation for AField<float,SIMD2>.
template<>
const std::string ASolver_BiCGStab_Cmplx<AField<float,SIMD2> >::class_name
                     = "ASolver_BiCGStab_Cmplx<AField<float,SIMD2> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_BiCGStab_Cmplx<AField_dev<float,SIMD2> >::register_factory();
}
#endif

template class ASolver_BiCGStab_Cmplx<AField<float,SIMD2> >;


//====================================================================
//============================================================END=====
