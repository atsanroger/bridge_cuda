/*!
      @file    asolver_CGNR.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib_alt/Solver/asolver_CGNR.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"

#include "lib_alt/Solver/asolver_CGNR-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,ACCEL>.
template<>
const std::string ASolver_CGNR<AField<double,ACCEL> >::class_name
                              = "ASolver_CGNR<Afield<double,ACCEL> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CGNR<AField<double,ACCEL> >::register_factory();
}
#endif

template class ASolver_CGNR<AField<double,ACCEL> >;

//====================================================================
// explicit instanciation for AField<float,ACCEL>.
template<>
const std::string ASolver_CGNR<AField<float,ACCEL> >::class_name
                              = "ASolver_CGNR<Afield<float,ACCEL> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CGNR<AField<float,ACCEL> >::register_factory();
}
#endif

template class ASolver_CGNR<AField<float,ACCEL> >;

//============================================================END=====
