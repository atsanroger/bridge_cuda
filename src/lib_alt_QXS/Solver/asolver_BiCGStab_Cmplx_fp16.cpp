/*!
        @file    asolver_BiCGStab_Cmplx_fp16.cpp
        @brief   BiCGStab Solver
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate: 2024-01-13 22:19:26 +0900 (2024年01月13日 (土)) $
        @version $LastChangedRevision: 2569 $
*/
//====================================================================

#include "lib_alt/Solver/asolver_BiCGStab_Cmplx.h"
#include "lib/ResourceManager/threadManager.h"
#include "lib_alt_QXS/Field/afield.h"

// function template files
#include "lib_alt_QXS/Field/afield-inc.h"

// template functions
#include "lib_alt/Solver/asolver_BiCGStab_Cmplx_fp16-tmpl.h"


//====================================================================
// explicit instanciation for AField<half>.
#ifdef USE_QXS_FP16

template<>
const std::string ASolver_BiCGStab_Cmplx<AField<half, QXS> >::class_name
  = "ASolver_BiCGStab_Cmplx<AField<half,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_BiCGStab_Cmplx<AField_dev<half, QXS> >::register_factory();
}
#endif

template class ASolver_BiCGStab_Cmplx<AField<half, QXS> >;

#endif
//============================================================END=====
