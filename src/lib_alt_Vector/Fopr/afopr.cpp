/*!
        @file    afopr.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
        @version $LastChangedRevision: 2543 $
*/

#include "lib/Fopr/afopr.h"
#include "lib_alt_Vector/Field/afield.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else

#include "lib_alt_Vector/Fopr/afopr_Wilson.h"
//#include "lib_alt_Vector/Fopr/afopr_Wilson_eo.h"
//#include "lib_alt_Vector/Fopr/afopr_Clover.h"
//#include "lib_alt_Vector/Fopr/afopr_Clover_eo.h"
#include "lib_alt_Vector/Fopr/afopr_Staggered.h"
//#include "lib_alt_Vector/Fopr/afopr_Staggered_eo.h"
#include "lib/Fopr/afopr_Domainwall.h"
//#include "lib_alt_Vector/Fopr/afopr_Domainwall_5din.h"
#include "lib/Fopr/afopr_Smeared.h"
//#include "lib/Fopr/afopr_Rational.h"


template<>
bool AFopr<AField<double,VECTOR> >::init_factory()
{
  typedef AField<double,VECTOR> AFIELD;
  bool result = true;
  result &= AFopr_Wilson<AFIELD>::register_factory();
  //  result &= AFopr_Wilson_eo<AFIELD>::register_factory();
  //  result &= AFopr_Clover<AFIELD>::register_factory();
  //  result &= AFopr_Clover_eo<AFIELD>::register_factory();
  result &= AFopr_Staggered<AFIELD>::register_factory();
  //  result &= AFopr_Staggered_eo<AFIELD>::register_factory();
  result &= AFopr_Smeared<AFIELD>::register_factory();
  //  result &= AFopr_Rational<AFIELD>::register_factory();
  result &= AFopr_Domainwall<AFIELD>::register_factory();
  //  result &= AFopr_Domainwall_5din<AFIELD>::register_factory();
  return result;
}

template<>
bool AFopr<AField<float,VECTOR> >::init_factory()
{
  typedef AField<float,VECTOR> AFIELD;
  bool result = true;
  result &= AFopr_Wilson<AFIELD>::register_factory();
  //  result &= AFopr_Wilson_eo<AFIELD>::register_factory();
  //  result &= AFopr_Clover<AFIELD>::register_factory();
  //  result &= AFopr_Clover_eo<AFIELD>::register_factory();
  result &= AFopr_Staggered<AFIELD>::register_factory();
  //  result &= AFopr_Staggered_eo<AFIELD>::register_factory();
  result &= AFopr_Smeared<AFIELD>::register_factory();
  //  result &= AFopr_Rational<AFIELD>::register_factory();
  result &= AFopr_Domainwall<AFIELD>::register_factory();
  //  result &= AFopr_Domainwall_5din<AFIELD>::register_factory();
  return result;
}

#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */


// explicit instanciation.
template class AFopr<AField<float,VECTOR> >;
template class AFopr<AField<double,VECTOR> >;

//============================================================END=====
