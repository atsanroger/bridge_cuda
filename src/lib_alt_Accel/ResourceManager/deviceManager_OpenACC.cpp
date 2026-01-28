/*!
        @file    deviceManager_OpenACC.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-08-20 14:25:12 #$
        @version $LastChangedRevision: 2535 $
*/
#ifdef USE_ACCEL_OPENACC

#include "ResourceManager/deviceManager_OpenACC.h"

// the following code only applies to NVIDIA
#ifdef _ACCEL
#include "accel.h"
#endif

#include "lib/Communicator/communicator.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

int DeviceManager_OpenACC::m_num_devices = 1;
int DeviceManager_OpenACC::m_device_id  = 0;
Bridge::VerboseLevel DeviceManager_OpenACC::m_vl = Bridge::CRUCIAL;
const std::string DeviceManager_OpenACC::class_name
                                           = "DeviceManager_OpenACC";
//====================================================================
void DeviceManager_OpenACC::init(int device_num_offset = 0)
{
  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: being initialized.\n", class_name.c_str());

  int process_id = Communicator::nodeid();
  vout.general(m_vl, "  process_id = %d\n", process_id);

#ifdef _ACCEL // the following code only applies to NVIDIA
  vout.general(m_vl, "  NVIDIA environment is used.\n");
  m_num_devices = acc_get_num_devices(acc_device_nvidia);
  m_device_id = process_id % m_num_devices + device_num_offset;
  acc_set_device_num(m_device_id, acc_device_nvidia);
#endif
  vout.general(m_vl, "  number of devices = %d\n", m_num_devices);
  // pprintf("  mydevice  = %d\n", mydevice);
  vout.general(m_vl, "  device_id = %d  is used on process = %d\n",
               m_device_id, process_id);

}

//====================================================================
void DeviceManager_OpenACC::finalize()
{
  vout.paranoiac(m_vl, "%s: finalization.\n", class_name.c_str());
}

#endif
//============================================================END=====
