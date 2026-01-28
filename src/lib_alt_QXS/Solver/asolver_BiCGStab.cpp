/*!
      @file    asolver_BiCGStab.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2569 $
*/

#include "lib_alt/Solver/asolver_BiCGStab.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"


//#include "lib_alt/Solver/asolver_BiCGStab-tmpl.h"
#include "lib_alt/Solver/asolver_BiCGStab_fp64-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,QXS>.
template<>
const std::string ASolver_BiCGStab<AField<double, QXS> >::class_name
  = "ASolver_BiCGStab<Afield<double,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_BiCGStab<AField_dev<double, QXS> >::register_factory();
}
#endif

template class ASolver_BiCGStab<AField<double, QXS> >;

//====================================================================
// explicit instanciation for AField<float,QXS>.
template<>
const std::string ASolver_BiCGStab<AField<float, QXS> >::class_name
  = "ASolver_BiCGStab<AField<float,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_BiCGStab<AField_dev<float, QXS> >::register_factory();
}
#endif

template class ASolver_BiCGStab<AField<float, QXS> >;


//============================================================END=====
