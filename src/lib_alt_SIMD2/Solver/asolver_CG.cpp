/*!
      @file    asolver_CG.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib_alt/Solver/asolver_CG.h"

#include "lib/ResourceManager/threadManager_OpenMP.h"

#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/afield-inc.h"

#include "lib_alt/Solver/asolver_CG-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,SIMD2>.
template<>
const std::string ASolver_CG<AField<double,SIMD2> >::class_name
                              = "ASolver_CG<Afield<double,SIMD2> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CG<AField<double,SIMD2> >::register_factory();
}
#endif

template class ASolver_CG<AField<double,SIMD2> >;

//====================================================================
// explicit instanciation for AField<float,SIMD2>.
template<>
const std::string ASolver_CG<AField<float,SIMD2> >::class_name
                              = "ASolver_CG<Afield<float,SIMD2> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CG<AField<float,SIMD2> >::register_factory();
}
#endif

template class ASolver_CG<AField<float,SIMD2> >;

//============================================================END=====
