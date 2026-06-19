/*!
      @file    afopr_Rational.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2527 $
*/

#include "lib/Fopr/afopr_Rational.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"

template<>
const std::string AFopr_Rational<AField<double, QXS> >::class_name
                            = "AFopr_Rational<AField<double,QXS> >";

template<>
const std::string AFopr_Rational<AField<float, QXS> >::class_name
                             = "AFopr_Rational<AField<float,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr<AField<double, QXS> >::Factory_params::Register(
    "Rational", create_object);

  bool init2 = AFopr<AField<float, QXS> >::Factory_params::Register(
    "Rational", create_object);
}
#endif

#include "lib/Fopr/afopr_Rational-tmpl.h"

// explicit instanciation.
template class AFopr_Rational<AField<double, QXS> >;
template class AFopr_Rational<AField<float, QXS> >;

//============================================================END=====
