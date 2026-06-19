/*!
        @file    init_alt_SIMD2.cpp
        @brief
        @author  Hideo Matsufuru  (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#include "lib_alt_SIMD2/init_alt_SIMD2.h"
#include "lib_alt_SIMD2/bridge_init_factory_alt_SIMD2.h"

bool init_alt_SIMD2()
{
  bool result = true;

  // Factory initialization
#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
  result &= bridge_init_factory_alt_SIMD2();
#endif
  
#ifdef DEBUG
  bridge_report_factory_alt_SIMD2();
#endif

#endif /* USE_FACTORY */

  return result;
}

bool fin_alt_SIMD2()
{
  return true;
}
