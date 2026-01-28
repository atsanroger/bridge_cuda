/*!
      @file    asolver_BiCGStab_fp16.cpp
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
#include "lib_alt/Solver/asolver_BiCGStab_fp16-tmpl.h"

//====================================================================
// explicit instanciation for AField<half,QXS>.
#ifdef USE_QXS_FP16

template<>
const std::string ASolver_BiCGStab<AField<half, QXS> >::class_name
  = "ASolver_BiCGStab<AField<half,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_BiCGStab<AField_dev<half, QXS> >::register_factory();
}
#endif

template class ASolver_BiCGStab<AField<half, QXS> >;

#endif
//============================================================END=====
