/*!                                                                             
        @file    afopr_Rational.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#include "lib/Fopr/afopr_Rational.h"
#include "lib/Fopr/afopr_Rational-tmpl.h"
#include "lib/Solver/ashiftsolver_CG.h"

#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/afield-inc.h"


template<>
const std::string AFopr_Rational<AField<double,SIMD2> >::class_name
= "AFopr_Rational<AField<double,SIMD2> >";

template<>
const std::string AFopr_Rational<AField<float,SIMD2> >::class_name
= "AFopr_Rational<AField<float,SIMD2> >";


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr<AField<double,SIMD2> >::Factory_params::Register(
                                          "Rational", create_object);

  //  bool init2 = AFopr<AField<float,SIMD2> >::Factory_params::Register(
  //                                          "Rational", create_object);
}
#endif


// explicit instanciation.
template class AFopr_Rational<AField<double,SIMD2> >;
//template class AFopr_Rational<AField<float,SIMD2> >;

//============================================================END=====
