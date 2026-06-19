/*!
      @file    asolver_CG.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib_alt/Solver/asolver_CG.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"

#include "lib_alt/Solver/asolver_CG-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,ACCEL>.
template<>
const std::string ASolver_CG<AField<double,ACCEL> >::class_name
                              = "ASolver_CG<Afield<double,ACCEL> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CG<AField<double,ACCEL> >::register_factory();
}
#endif

template class ASolver_CG<AField<double,ACCEL> >;

//====================================================================
// explicit instanciation for AField<float,ACCEL>.
template<>
const std::string ASolver_CG<AField<float,ACCEL> >::class_name
                              = "ASolver_CG<Afield<float,ACCEL> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CG<AField<float,ACCEL> >::register_factory();
}
#endif

template class ASolver_CG<AField<float,ACCEL> >;

//============================================================END=====
