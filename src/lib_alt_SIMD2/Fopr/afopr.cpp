/*!                                                                             
        @file    afopr.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#include "lib/Fopr/afopr.h"
#include "lib_alt_SIMD2/Field/afield.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else

#include "lib_alt_SIMD2/Fopr/afopr_Wilson.h"
#include "lib_alt_SIMD2/Fopr/afopr_Wilson_eo.h"
#include "lib_alt_SIMD2/Fopr/afopr_Clover.h"
#include "lib_alt_SIMD2/Fopr/afopr_Clover_eo.h"
#include "lib_alt_SIMD2/Fopr/afopr_Staggered.h"
//#include "lib_alt_SIMD2/Fopr/afopr_Staggered_eo.h"
#include "lib_alt/Fopr/afopr_Domainwall_General.h"
#include "lib_alt/Fopr/afopr_Domainwall_General_eo.h"
#include "lib_alt_SIMD2/Fopr/afopr_Domainwall_General_5din.h"
#include "lib_alt_SIMD2/Fopr/afopr_Domainwall_General_5din_eo.h"
#include "lib/Fopr/afopr_Smeared.h"
#include "lib/Fopr/afopr_Rational.h"


template<>
bool AFopr<AField<double,SIMD2> >::init_factory()
{
  typedef AField<double,SIMD2> AFIELD;
  bool result = true;
  result &= AFopr_Wilson<AFIELD>::register_factory();
  result &= AFopr_Wilson_eo<AFIELD>::register_factory();
  result &= AFopr_Clover<AFIELD>::register_factory();
  result &= AFopr_Clover_eo<AFIELD>::register_factory();
  result &= AFopr_Staggered<AFIELD>::register_factory();
  //  result &= AFopr_Staggered_eo<AFIELD>::register_factory();
  result &= AFopr_Smeared<AFIELD>::register_factory();
  result &= AFopr_Rational<AFIELD>::register_factory();
  result &= AFopr_Domainwall_General<AFIELD>::register_factory();
  result &= AFopr_Domainwall_General_eo<AFIELD>::register_factory();
  result &= AFopr_Domainwall_General_5din<AFIELD>::register_factory();
  result &= AFopr_Domainwall_General_5din_eo<AFIELD>::register_factory();
  return result;
}

template<>
bool AFopr<AField<float,SIMD2> >::init_factory()
{
  typedef AField<float,SIMD2> AFIELD;
  bool result = true;
  result &= AFopr_Wilson<AFIELD>::register_factory();
  result &= AFopr_Wilson_eo<AFIELD>::register_factory();
  result &= AFopr_Clover<AFIELD>::register_factory();
  result &= AFopr_Clover_eo<AFIELD>::register_factory();
  result &= AFopr_Staggered<AFIELD>::register_factory();
  //  result &= AFopr_Staggered_eo<AFIELD>::register_factory();
  result &= AFopr_Smeared<AFIELD>::register_factory();
  //  result &= AFopr_Rational<AFIELD>::register_factory();
  result &= AFopr_Domainwall_General<AFIELD>::register_factory();
  result &= AFopr_Domainwall_General_eo<AFIELD>::register_factory();
  result &= AFopr_Domainwall_General_5din<AFIELD>::register_factory();
  result &= AFopr_Domainwall_General_5din_eo<AFIELD>::register_factory();
  return result;
}

#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */


template<>
const std::string AFopr<AField<double,SIMD2> >::
    class_name = "AFopr<AField<double,SIMD2> >";

template<>
const std::string AFopr<AField<float,SIMD2> >::
    class_name = "AFopr<AField<float,SIMD2> >";


// explicit instanciation.
template class AFopr<AField<float,SIMD2> >;
template class AFopr<AField<double,SIMD2> >;

//============================================================END=====
