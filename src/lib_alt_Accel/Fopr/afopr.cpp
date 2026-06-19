/*!                                                                            
        @file    afopr.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2024-03-26 15:47:37 #$
        @version $LastChangedRevision: 2590 $
*/

#include "lib/Fopr/afopr.h"
#include "lib_alt_Accel/Field/afield.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else

#include "lib_alt_Accel/Fopr/afopr_Wilson.h"
#include "lib_alt_Accel/Fopr/afopr_Wilson_eo.h"
#include "lib_alt_Accel/Fopr/afopr_Clover.h"
#include "lib_alt_Accel/Fopr/afopr_Wilson_dd.h"
#include "lib_alt_Accel/Fopr/afopr_Clover_dd.h"
#include "lib_alt_Accel/Fopr/afopr_Clover_eo.h"
#include "lib_alt_Accel/Fopr/afopr_Staggered.h"
#include "lib_alt_Accel/Fopr/afopr_Staggered_eo.h"
#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din.h"
#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din_dd.h"
#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din_eo.h"
#include "lib_alt_Accel/Fopr/afopr_Domainwall_PVprec.h"
#include "lib_alt_Accel/Fopr/afopr_Staggered_CLE.h"
#include "lib_alt_Accel/Fopr/afopr_Wilson_CLE.h"
#include "lib/Fopr/afopr_Domainwall.h"
#include "lib/Fopr/afopr_Domainwall_eo.h"
//#include "lib/Fopr/afopr_OptimalDomainwall.h"
#include "lib/Fopr/afopr_Smeared.h"
#include "lib/Fopr/afopr_Rational.h"


template<>
bool AFopr<AField<double,ACCEL> >::init_factory()
{
  typedef AField<double,ACCEL> AFIELD;
  bool result = true;
  result &= AFopr_Wilson<AFIELD>::register_factory();
  result &= AFopr_Clover<AFIELD>::register_factory();
  result &= AFopr_Wilson_dd<AFIELD>::register_factory();
  result &= AFopr_Clover_dd<AFIELD>::register_factory();
  result &= AFopr_Domainwall_5din<AFIELD>::register_factory();
  result &= AFopr_Domainwall_5din_dd<AFIELD>::register_factory();
  result &= AFopr_Domainwall_5din_eo<AFIELD>::register_factory();
  result &= AFopr_Domainwall_PVprec<AFIELD>::register_factory();
  result &= AFopr_Domainwall<AFIELD>::register_factory();
  result &= AFopr_Domainwall_eo<AFIELD>::register_factory();
#ifdef USE_ACCEL_OPENACC
  result &= AFopr_Wilson_eo<AFIELD>::register_factory();
  result &= AFopr_Clover_eo<AFIELD>::register_factory();
  result &= AFopr_Staggered<AFIELD>::register_factory();
  result &= AFopr_Staggered_eo<AFIELD>::register_factory();
  result &= AFopr_Smeared<AFIELD>::register_factory();
  result &= AFopr_Rational<AFIELD>::register_factory();
  //  result &= AFopr_OptimalDomainwall<AFIELD>::register_factory();
  result &= AFopr_Staggered_CLE<AFIELD>::register_factory();
  result &= AFopr_Wilson_CLE<AFIELD>::register_factory();
#endif
  return result;
}

template<>
bool AFopr<AField<float,ACCEL> >::init_factory()
{
  typedef AField<float,ACCEL> AFIELD;
  bool result = true;
  result &= AFopr_Wilson<AFIELD>::register_factory();
  result &= AFopr_Clover<AFIELD>::register_factory();
  result &= AFopr_Wilson_dd<AFIELD>::register_factory();
  result &= AFopr_Clover_dd<AFIELD>::register_factory();
  result &= AFopr_Domainwall_5din<AFIELD>::register_factory();
  result &= AFopr_Domainwall_5din_dd<AFIELD>::register_factory();
  result &= AFopr_Domainwall_5din_eo<AFIELD>::register_factory();
  result &= AFopr_Domainwall_PVprec<AFIELD>::register_factory();
  result &= AFopr_Domainwall<AFIELD>::register_factory();
  result &= AFopr_Domainwall_eo<AFIELD>::register_factory();
#ifdef USE_ACCEL_OPENACC
  result &= AFopr_Wilson_eo<AFIELD>::register_factory();
  result &= AFopr_Clover_eo<AFIELD>::register_factory();
  result &= AFopr_Staggered<AFIELD>::register_factory();
  result &= AFopr_Staggered_eo<AFIELD>::register_factory();
  result &= AFopr_Smeared<AFIELD>::register_factory();
  //  result &= AFopr_Rational<AFIELD>::register_factory();
  //  result &= AFopr_OptimalDomainwall<AFIELD>::register_factory();
#endif
  return result;
}

#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */


// explicit instanciation.
template class AFopr<AField<double,ACCEL> >;
template class AFopr<AField<float,ACCEL> >;

//============================================================END=====
