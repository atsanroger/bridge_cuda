/*!
        @file    deviceManager_Vector.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
        @version $LastChangedRevision: 2543 $
*/

#include "lib_alt_Vector/ResourceManager/deviceManager_Vector.h"

#include "lib/Communicator/communicator.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

int DeviceManager_Vector::m_num_devices = 1;
int DeviceManager_Vector::m_device_id  = 0;
Bridge::VerboseLevel DeviceManager_Vector::m_vl = Bridge::CRUCIAL;
const std::string DeviceManager_Vector::class_name
                                           = "DeviceManager_Vector";
//====================================================================
void DeviceManager_Vector::init()
{
  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: being initialized.\n", class_name.c_str());

  int process_id = Communicator::nodeid();
  vout.general(m_vl, "  process_id = %d\n", process_id);

}

//====================================================================
void DeviceManager_Vector::finalize()
{
  vout.paranoiac(m_vl, "%s: finalization.\n", class_name.c_str());
}

//====================================================================
//============================================================END=====
