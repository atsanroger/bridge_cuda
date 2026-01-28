/*!
        @file    deviceManager_CUDA.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-09-12 13:12:27 #$
        @version $LastChangedRevision: 2539 $
*/

#ifdef USE_ACCEL_CUDA

#include "lib_alt_Accel/ResourceManager/deviceManager_CUDA.h"

// the following code only applies to NVIDIA
#ifdef _ACCEL
#include "accel.h"
//#include "cuda_runtime.h"
#endif

#include "lib/Communicator/communicator.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

int DeviceManager_CUDA::m_num_devices = 1;
int DeviceManager_CUDA::m_device_id  = 0;

Bridge::VerboseLevel DeviceManager_CUDA::m_vl = Bridge::CRUCIAL;

const std::string DeviceManager_CUDA::class_name
                                           = "DeviceManager_CUDA";

//====================================================================
void DeviceManager_CUDA::init(int device_num_offset = 0)
{
  m_vl = CommonParameters::Vlevel();
  vout.general(m_vl, "%s: being initialized.\n", class_name.c_str());

  int process_id = Communicator::nodeid();
  vout.general(m_vl, "  process_id = %d\n", process_id);

#ifdef _ACCEL
  vout.general(m_vl, "  NVIDIA environment is used.\n");
  cudaGetDeviceCount(&m_num_devices);

  m_device_id = process_id % m_num_devices + device_num_offset;
  cudaSetDevice(m_device_id);

  for (int dev = 0; dev < m_num_devices; ++dev)
  {
    if (dev == m_device_id) {
      continue;
    }
    int can_access = 0;
    cudaDeviceCanAccessPeer(&can_access, m_device_id, dev);
    if (can_access == 1) {
      cudaDeviceEnablePeerAccess(dev, 0);
    }
  }
#endif

  vout.general(m_vl, "  number of devices = %d\n", m_num_devices);
  vout.general(m_vl, "  device_id = %d  is used on process = %d\n",
               m_device_id, process_id);
}


//====================================================================
void DeviceManager_CUDA::finalize()
{
  vout.paranoiac(m_vl, "%s: finalization.\n", class_name.c_str());
}


#endif
//============================================================END=====
