/*!
      @file    afopr_Precond_dd_float.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2604 $
*/


#include "lib/ResourceManager/threadManager.h"
#include "complexTraits.h"

#include "afopr_Precond_dd.h"

#ifdef USE_ALT_QXS

#include "lib_alt_QXS/Field/afield.h"

typedef AField<float, QXS> AField_f;

template<>
const std::string AFopr_Precond_dd<AField_f >::class_name
                     = "AFopr_Precond_dd<AField<float,QXS> >";

#endif


#ifdef USE_ALT_ACCEL

#include "lib_alt_Accel/Field/afield.h"

typedef AField<float, ACCEL> AField_f;

template<>
const std::string AFopr_Precond_dd<AField_f>::class_name
                     = "AFopr_Precond_dd<AField<float,ACCEL> >";

#endif


#ifdef USE_FACTORY_AUTOREGISTER
/*
namespace {
  bool init2 = AFopr<AField_f >::Factory_params::Register(
                          "Precond_dd", create_object_with_params);
}
*/
#endif


#include "afopr_Precond_dd-tmpl.h"

// explicit instanciation

template class AFopr_Precond_dd<AField_f>;


//============================================================END=====
