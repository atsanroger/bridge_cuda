/*!
      @file    afopr_Smeared.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib/Fopr/afopr_Smeared.h"
#include "lib/Fopr/afopr_Smeared-tmpl.h"

#include "lib_alt_SIMD2/Field/afield.h"


template<>
const std::string AFopr_Smeared<AField<double,SIMD2> >::class_name
                                = "AFopr_Smeared<AField<double,SIMD2> >";

template<>
const std::string AFopr_Smeared<AField<float,SIMD2> >::class_name
                                = "AFopr_Smeared<AField<float,SIMD2> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr<AField<double,SIMD2> >::Factory_params::Register(
                                         "Smeared", create_object);

  bool init2 = AFopr<AField<float,SIMD2> >::Factory_params::Register(
                                         "Smeared", create_object);
}
#endif

// explicit instanciation.
template class AFopr_Smeared<AField<float,SIMD2> >;
template class AFopr_Smeared<AField<double,SIMD2> >;

//============================================================END=====
