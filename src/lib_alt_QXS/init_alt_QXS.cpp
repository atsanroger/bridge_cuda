/*!
        @file    init_alt_QXS.cpp
        @brief
        @author  Hideo Matsufuru  (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
        @version $LastChangedRevision: 2499 $
*/

#include "lib_alt_QXS/init_alt_QXS.h"
#include "lib_alt_QXS/bridge_init_factory_alt_QXS.h"

bool init_alt_QXS()
{
  bool result = true;

  // Factory initialization
#ifdef USE_FACTORY
#ifdef USE_FACTORY_AUTOREGISTER
#else
  result &= bridge_init_factory_alt_QXS();
#endif

#ifdef DEBUG
  bridge_report_factory_alt_QXS();
#endif
#endif /* USE_FACTORY */

  return result;
}


bool fin_alt_QXS()
{
  bool result = true;

  // Factory finalization
#ifdef USE_FACTORY
  //result &= bridge_fin_factory_alt_QXS();
#endif
  return true;
}
