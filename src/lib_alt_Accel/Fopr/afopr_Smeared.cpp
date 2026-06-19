/*!
      @file    afopr_Smeared.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib/Fopr/afopr_Smeared.h"
#include "lib/Fopr/afopr_Smeared-tmpl.h"

#include "lib_alt_Accel/Field/afield.h"


template<>
const std::string AFopr_Smeared<AField<double,ACCEL> >::class_name
                                = "AFopr_Smeared<AField<double,ACCEL> >";

template<>
const std::string AFopr_Smeared<AField<float,ACCEL> >::class_name
                                = "AFopr_Smeared<AField<float,ACCEL> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr<AField<double,ACCEL> >::Factory_params::Register(
                                         "Smeared", create_object);

  bool init2 = AFopr<AField<float,ACCEL> >::Factory_params::Register(
                                         "Smeared", create_object);
}
#endif

// explicit instanciation.
template class AFopr_Smeared<AField<float,ACCEL> >;
template class AFopr_Smeared<AField<double,ACCEL> >;

//============================================================END=====
