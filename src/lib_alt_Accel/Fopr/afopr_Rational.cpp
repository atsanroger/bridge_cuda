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

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"


template<>
const std::string AFopr_Rational<AField<double,ACCEL> >::class_name
= "AFopr_Rational<AField<double,ACCEL> >";

template<>
const std::string AFopr_Rational<AField<float,ACCEL> >::class_name
= "AFopr_Rational<AField<float,ACCEL> >";


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr<AField<double,ACCEL> >::Factory_params::Register(
                                          "Rational", create_object);

  //  bool init2 = AFopr<AField<float,ACCEL> >::Factory_params::Register(
  //                                          "Rational", create_object);
}
#endif


// explicit instanciation.
template class AFopr_Rational<AField<double,ACCEL> >;
//template class AFopr_Rational<AField<float,ACCEL> >;

//============================================================END=====
