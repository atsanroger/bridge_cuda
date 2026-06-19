/*!
      @file    afopr_Smeared.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#include "lib/Fopr/afopr_Smeared.h"
#include "lib/Fopr/afopr_Smeared-tmpl.h"

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/Field/afield.h"


template<>
const std::string AFopr_Smeared<AField<double, QXS> >::class_name
  = "AFopr_Smeared<AField<double,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr<AField<double, QXS> >::Factory_params::Register(
    "Smeared", create_object);
}
#endif

// explicit instanciation.
template class AFopr_Smeared<AField<double, QXS> >;

//===================================================================
template<>
const std::string AFopr_Smeared<AField<float, QXS> >::class_name
  = "AFopr_Smeared<AField<float,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AFopr<AField<float, QXS> >::Factory_params::Register(
    "Smeared", create_object);
}
#endif

// explicit instanciation.
template class AFopr_Smeared<AField<float, QXS> >;


//===================================================================
#ifdef USE_QXS_FP16

template<>
const std::string AFopr_Smeared<AField<half, QXS> >::class_name
  = "AFopr_Smeared<AField<half,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AFopr<AField<half, QXS> >::Factory_params::Register(
    "Smeared", create_object);
}
#endif

// explicit instanciation.
template class AFopr_Smeared<AField<half, QXS> >;

#endif
//============================================================END=====
