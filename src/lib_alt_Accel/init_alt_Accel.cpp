/*!
        @file    init_alt_Accel.cpp
        @brief
        @author  Hideo Matsufuru  (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2024-10-27 21:04:29 #$
        @version $LastChangedRevision: 2610 $
*/

#include "lib_alt_Accel/init_alt_Accel.h"
#include "lib_alt_Accel/bridge_init_factory_alt_Accel.h"

#ifdef USE_ACCEL_OPENACC
#include "lib_alt_Accel/ResourceManager/deviceManager_OpenACC.h"
#endif

#ifdef USE_ACCEL_CUDA
#include "lib_alt_Accel/ResourceManager/deviceManager_CUDA.h"
#endif

#ifdef USE_ACCEL_NVSHMEM
#include "lib_alt_Accel/ResourceManager/deviceManager_NVSHMEM.h"
#endif

//====================================================================
bool init_alt_Accel(Parameters& params)
{
  vout.general("Alt_Accel code is being setup.\n");

  bool result = true;

  int device_num_offset = 0;
  bool set_offset = params.is_set("device_number_offset");
  if(set_offset) device_num_offset
                  = params.get_int("device_number_offset");

#ifdef USE_ACCEL_OPENACC
  DeviceManager_OpenACC::init(device_num_offset);
#endif

#ifdef USE_ACCEL_CUDA
  DeviceManager_CUDA::init(device_num_offset);
#endif

#ifdef USE_ACCEL_NVSHMEM
  DeviceManager_NVSHMEM::init(device_num_offset);
#endif

// Factory initialization
#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
  result &= bridge_init_factory_alt_Accel();
#endif
  
bridge_report_factory_alt_Accel();

#endif /* USE_FACTORY */

  return result;

}

//====================================================================
bool fin_alt_Accel()
{
  return true;
}

//============================================================END=====
