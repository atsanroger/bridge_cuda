/*!
      @file    asolver_CG_fp16.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2562 $
*/

#include "lib_alt/Solver/asolver_CG.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"

#include "lib_alt/Solver/asolver_CG_fp16-tmpl.h"
//#include "lib_alt/Solver/asolver_CG-tmpl.h"

//====================================================================
// explicit instanciation for AField<half,QXS>.
#ifdef USE_QXS_FP16

template<>
const std::string ASolver_CG<AField<half, QXS> >::class_name
  = "ASolver_CG<Afield<half,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CG<AField<half, QXS> >::register_factory();
}
#endif

template class ASolver_CG<AField<half, QXS> >;

#endif
//============================================================END=====
