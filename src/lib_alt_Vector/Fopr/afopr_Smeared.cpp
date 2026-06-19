/*!
      @file    afopr_Smeared.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

#include "lib/Fopr/afopr_Smeared.h"
#include "lib/Fopr/afopr_Smeared-tmpl.h"

#include "lib_alt_Vector/Field/afield.h"


template<>
const std::string AFopr_Smeared<AField<double,VECTOR> >::class_name
                                = "AFopr_Smeared<AField<double,VECTOR> >";

template<>
const std::string AFopr_Smeared<AField<float,VECTOR> >::class_name
                                = "AFopr_Smeared<AField<float,VECTOR> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr<AField<double,VECTOR> >::Factory_params::Register(
                                         "Smeared", create_object);

  bool init2 = AFopr<AField<float,VECTOR> >::Factory_params::Register(
                                         "Smeared", create_object);
}
#endif

// explicit instanciation.
template class AFopr_Smeared<AField<float,VECTOR> >;
template class AFopr_Smeared<AField<double,VECTOR> >;

//============================================================END=====
