/*!
        @file    fopr_Rational.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2021-06-15 22:41:26 #$

        @version $LastChangedRevision: 2271 $
*/

#include "Fopr/fopr_Rational.h"
#include "Fopr/afopr_Rational-tmpl.h"
#include "Field/field.h"


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Rational::register_factory();
}
#endif

template<>
const std::string AFopr_Rational<Field>::class_name = "Fopr_Rational";

template class AFopr_Rational<Field>;

//====================================================================
//============================================================END=====
