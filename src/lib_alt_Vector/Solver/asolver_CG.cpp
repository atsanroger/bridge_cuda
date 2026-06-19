/*!
      @file    asolver_CG.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

#include "lib_alt/Solver/asolver_CG.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Vector/Field/afield.h"
#include "lib_alt_Vector/Field/afield-inc.h"

#include "lib_alt/Solver/asolver_CG-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,VECTOR>.
template<>
const std::string ASolver_CG<AField<double,VECTOR> >::class_name
                              = "ASolver_CG<Afield<double,VECTOR> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CG<AField<double,VECTOR> >::register_factory();
}
#endif

template class ASolver_CG<AField<double,VECTOR> >;

//====================================================================
// explicit instanciation for AField<float,VECTOR>.
template<>
const std::string ASolver_CG<AField<float,VECTOR> >::class_name
                              = "ASolver_CG<Afield<float,VECTOR> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CG<AField<float,VECTOR> >::register_factory();
}
#endif

template class ASolver_CG<AField<float,VECTOR> >;

//============================================================END=====
