/*!
        @file    init_alt_Vector.cpp
        @brief
        @author  Hideo Matsufuru  (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
        @version $LastChangedRevision: 2543 $
*/

#include "lib_alt_Vector/init_alt_Vector.h"
#include "lib_alt_Vector/ResourceManager/deviceManager_Vector.h"
#include "lib_alt_Vector/bridge_init_factory_alt_Vector.h"

bool init_alt_Vector()
{
  bool result = true;

  DeviceManager_Vector::init();

  // Factory initialization
#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
  result &= bridge_init_factory_alt_Vector();
#endif
  
#ifdef DEBUG
  bridge_report_factory_alt_Vector();
#endif

#endif /* USE_FACTORY */

  return result;
}

bool fin_alt_Vector()
{
  return true;
}
