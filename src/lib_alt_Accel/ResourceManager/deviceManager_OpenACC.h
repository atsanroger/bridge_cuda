/*!
        @file    deviceManager_OpenACC.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#ifndef DEVICEMANAGER_OPENACC_INCLUDED
#define DEVICEMANAGER_OPENACC_INCLUDED

#include <string>

#include "lib/Parameters/commonParameters.h"
#include "lib/IO/bridgeIO.h"

//! Device manger for OpenACC.

/*!
   Device manger for OpenACC.
   This class is static.
                                    [06 Jun 2019 H.Matsufuru]
 */
class DeviceManager_OpenACC {

 public:
  static const std::string class_name;

 private:
  static int m_num_devices;          //!< number of devices.
  static int m_device_id;            //!< device id used in this process.
  static Bridge::VerboseLevel m_vl;  //!< verbose level.

 public:
  //! setup: called in main only once.
  static void init(int device_num_offset);

  //! finalization.
  static void finalize();

 private:
  // non-copyable
  DeviceManager_OpenACC(const DeviceManager_OpenACC&);
  DeviceManager_OpenACC& operator=(const DeviceManager_OpenACC&);

};
#endif
