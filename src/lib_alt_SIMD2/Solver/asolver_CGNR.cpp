/*!
      @file    asolver_CGNR.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib_alt/Solver/asolver_CGNR.h"

#include "lib/ResourceManager/threadManager_OpenMP.h"

#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/afield-inc.h"

#include "lib_alt/Solver/asolver_CGNR-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,SIMD2>.
template<>
const std::string ASolver_CGNR<AField<double,SIMD2> >::class_name
                              = "ASolver_CGNR<Afield<double,SIMD2> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CGNR<AField<double,SIMD2> >::register_factory();
}
#endif

template class ASolver_CGNR<AField<double,SIMD2> >;

//====================================================================
// explicit instanciation for AField<float,SIMD2>.
template<>
const std::string ASolver_CGNR<AField<float,SIMD2> >::class_name
                              = "ASolver_CGNR<Afield<float,SIMD2> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CGNR<AField<float,SIMD2> >::register_factory();
}
#endif

template class ASolver_CGNR<AField<float,SIMD2> >;

//============================================================END=====
