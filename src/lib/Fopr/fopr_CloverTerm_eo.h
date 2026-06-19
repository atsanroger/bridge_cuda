/*!
        @file    fopr_CloverTerm_eo.h

        @brief

        @author  UEDA, Satoru (sueda)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef FOPR_CLOVERTERM_EO_INCLUDED
#define FOPR_CLOVERTERM_EO_INCLUDED

// implementations

// original unoptimised version
#ifdef USE_ORG
#include "Org/fopr_CloverTerm_eo_impl.h"
#endif

// improved version
#ifdef USE_IMP
#include "Imp/fopr_CloverTerm_eo_impl.h"
#endif

#if defined(USE_IMP)
typedef Imp::Fopr_CloverTerm_eo   Fopr_CloverTerm_eo;
#else
typedef Org::Fopr_CloverTerm_eo   Fopr_CloverTerm_eo;
#endif

#endif
