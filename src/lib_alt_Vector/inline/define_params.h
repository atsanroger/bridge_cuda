/*!
        @file    define_params_SU3.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-11-13 13:02:59 #$
        @version $LastChangedRevision: 2553 $
*/

#ifndef DEFINE_PARAMS_INCLUDED
#define DEFINE_PARAMS_INCLUDED

#include "bridge_defs.h"

#if defined USE_GROUP_SU3
#include "lib_alt_Vector/inline/define_params_SU3.h"
#endif

namespace {

  inline int ceil_nwp(const int nst)
  {
    int size = nst/NWP;
    if( (nst % NWP) > 0) ++size;
    return size * NWP;
  }

}


#endif
//============================================================END=====
