/*!
      @file    asolver_BiCGStab_Cmplx.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

#include "lib_alt/Solver/asolver_BiCGStab_Cmplx.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Vector/Field/afield.h"
#include "lib_alt_Vector/Field/afield-inc.h"

#include "lib_alt/Solver/asolver_BiCGStab_Cmplx-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,VECTOR>.
template<>
const std::string ASolver_BiCGStab_Cmplx<AField<double,VECTOR> >::class_name
                     = "ASolver_BiCGStab_Cmplx<Afield<double,VECTOR> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_BiCGStab_Cmplx<AField_dev<double,VECTOR> >::register_factory();
}
#endif

template class ASolver_BiCGStab_Cmplx<AField<double,VECTOR> >;

//====================================================================
// explicit instanciation for AField<float,VECTOR>.
template<>
const std::string ASolver_BiCGStab_Cmplx<AField<float,VECTOR> >::class_name
                     = "ASolver_BiCGStab_Cmplx<AField<float,VECTOR> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_BiCGStab_Cmplx<AField_dev<float,VECTOR> >::register_factory();
}
#endif

template class ASolver_BiCGStab_Cmplx<AField<float,VECTOR> >;


//====================================================================
//============================================================END=====
