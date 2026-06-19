/*!
        @file    fopr.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-11-13 13:02:59 #$

        @version $LastChangedRevision: 2553 $
*/

#include "Fopr/fopr.h"
#include "Fopr/fopr_eo.h"
#include "Field/field.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "fopr_Wilson.h"
#include "Org/fopr_Wilson_impl.h"
#include "Imp/fopr_Wilson_impl.h"
#include "fopr_Wilson_eo.h"
#include "Org/fopr_Wilson_eo_impl.h"
#include "Imp/fopr_Wilson_eo_impl.h"
#include "fopr_Clover.h"
#include "fopr_Clover_eo.h"
#include "fopr_Domainwall.h"
#include "fopr_Overlap.h"
#include "fopr_Sign.h"
#include "fopr_Rational.h"
#include "fopr_Smeared.h"
#include "fopr_Smeared_alt.h"
#include "fopr_Smeared_eo.h"
#include "fopr_Chebyshev.h"
#include "fopr_Wilson_TwistedMass.h"
#include "fopr_WilsonGeneral.h"
#include "Org/fopr_WilsonGeneral_impl.h"
#include "Imp/fopr_WilsonGeneral_impl.h"
#include "fopr_CloverGeneral.h"
#include "fopr_Wilson_Chemical.h"
#include "fopr_Clover_Chemical.h"
#include "fopr_Staggered.h"
#include "fopr_Staggered_eo.h"
#include "fopr_Wilson_SF.h"
#include "fopr_Clover_SF.h"
#include "fopr_Rational_SF.h"
#include "fopr_NonRelativistic.h"
#include "fopr_CRS.h"

template<>
bool Fopr::init_factory()
{
  bool result = true;
  result &= Selector_Fopr_Wilson::register_factory();
  result &= Org::Fopr_Wilson::register_factory();
  result &= Imp::Fopr_Wilson::register_factory();
  result &= Selector_Fopr_Wilson_eo::register_factory();
  result &= Org::Fopr_Wilson_eo::register_factory();
  result &= Imp::Fopr_Wilson_eo::register_factory();
  result &= Fopr_Clover::register_factory();
  result &= Fopr_Clover_eo::register_factory();
  result &= Fopr_Domainwall::register_factory();
  result &= Fopr_Overlap::register_factory();
  result &= Fopr_Sign::register_factory();
  result &= Fopr_Rational::register_factory();
  result &= Fopr_Smeared::register_factory();
#ifdef USE_ALT_CODE
  result &= Fopr_Smeared_alt::register_factory();
#endif
  result &= Fopr_Smeared_eo::register_factory();
  result &= Fopr_Chebyshev::register_factory();
  result &= Fopr_Wilson_TwistedMass::register_factory();
  result &= Selector_Fopr_WilsonGeneral::register_factory();
  result &= Org::Fopr_WilsonGeneral::register_factory();
  result &= Imp::Fopr_WilsonGeneral::register_factory();
  result &= Fopr_CloverGeneral::register_factory();
  result &= Fopr_Wilson_Chemical::register_factory();
  result &= Fopr_Clover_Chemical::register_factory();
  result &= Fopr_Staggered::register_factory();
  result &= Fopr_Staggered_eo::register_factory();
  result &= Fopr_Wilson_SF::register_factory();
  result &= Fopr_Clover_SF::register_factory();
  result &= Fopr_Rational_SF::register_factory();
  result &= Fopr_NonRelativistic::register_factory();
  result &= Fopr_CRS::register_factory();

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */

// explicit instanciation.
template<>
const std::string AFopr<Field>::class_name = "Fopr";
template class AFopr<Field>;

template<>
const std::string AFopr_eo<Field>::class_name = "Fopr_eo";
template class AFopr_eo<Field>;
