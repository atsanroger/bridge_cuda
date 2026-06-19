/*!
        @file    define_params.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2022-12-16 15:57:38 #$
        @version $LastChangedRevision: 2422 $
*/

#ifndef DEFINE_PARAMS_INCLUDED
#define DEFINE_PARAMS_INCLUDED

// QXS is available only for the SU(3) case.
#if defined USE_GROUP_SU3
#include "lib_alt_QXS/inline/define_params_SU3.h"
#endif

//#define USE_QXSLIB  // now this macro is defined in Makefile

//#define USE_QXSLIB  // comment out if QXS library is not used.

//#define USE_BENCHMARK  // if BridgeQXS library separately complied

#endif
//============================================================END=====
