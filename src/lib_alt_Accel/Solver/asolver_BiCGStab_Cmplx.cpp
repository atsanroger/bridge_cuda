/*!
      @file    asolver_BiCGStab_Cmplx.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib_alt/Solver/asolver_BiCGStab_Cmplx.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"

#include "lib_alt/Solver/asolver_BiCGStab_Cmplx-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,ACCEL>.
template<>
const std::string ASolver_BiCGStab_Cmplx<AField<double,ACCEL> >::class_name
                     = "ASolver_BiCGStab_Cmplx<Afield<double,ACCEL> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_BiCGStab_Cmplx<AField_dev<double,ACCEL> >::register_factory();
}
#endif

template class ASolver_BiCGStab_Cmplx<AField<double,ACCEL> >;

//====================================================================
// explicit instanciation for AField<float,ACCEL>.
template<>
const std::string ASolver_BiCGStab_Cmplx<AField<float,ACCEL> >::class_name
                     = "ASolver_BiCGStab_Cmplx<AField<float,ACCEL> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_BiCGStab_Cmplx<AField_dev<float,ACCEL> >::register_factory();
}
#endif

template class ASolver_BiCGStab_Cmplx<AField<float,ACCEL> >;

//============================================================END=====
