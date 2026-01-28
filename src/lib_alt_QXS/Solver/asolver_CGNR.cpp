/*!
      @file    asolver_CGNR.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2562 $
*/

#include "lib_alt/Solver/asolver_CGNR.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"

#include "lib_alt/Solver/asolver_CGNR-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,QXS>.
template<>
const std::string ASolver_CGNR<AField<double, QXS> >::class_name
  = "ASolver_CGNR<Afield<double,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CGNR<AField<double, QXS> >::register_factory();
}
#endif

template class ASolver_CGNR<AField<double, QXS> >;

//====================================================================
// explicit instanciation for AField<float,QXS>.
template<>
const std::string ASolver_CGNR<AField<float, QXS> >::class_name
  = "ASolver_CGNR<Afield<float,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CGNR<AField<float, QXS> >::register_factory();
}
#endif

template class ASolver_CGNR<AField<float, QXS> >;

//====================================================================
// explicit instanciation for AField<half,QXS>.
#ifdef USE_QXS_FP16

template<>
const std::string ASolver_CGNR<AField<half, QXS> >::class_name
  = "ASolver_CGNR<Afield<half,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CGNR<AField<half, QXS> >::register_factory();
}
#endif

template class ASolver_CGNR<AField<half, QXS> >;

#endif
//============================================================END=====
