/*!
        @file    fopr_WilsonGeneral.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef FOPR_WILSONGENERAL_INCLUDED
#define FOPR_WILSONGENERAL_INCLUDED

// implementations

// original unoptimised version
#ifdef USE_ORG
#include "Org/fopr_WilsonGeneral_impl.h"
#endif

// improved version
#ifdef USE_IMP
#include "Imp/fopr_WilsonGeneral_impl.h"
#endif

#if defined(USE_IMP)
typedef Imp::Fopr_WilsonGeneral   Fopr_WilsonGeneral;
#else
typedef Org::Fopr_WilsonGeneral   Fopr_WilsonGeneral;
#endif

#ifdef USE_FACTORY
namespace Selector_Fopr_WilsonGeneral
{
  bool register_factory();
}
#endif

#endif /* FOPR_WILSONGENERAL_INCLUDED */
