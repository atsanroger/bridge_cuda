/*!
        @file    deviceManager_CUDA.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-08-20 14:25:12 #$
        @version $LastChangedRevision: 2535 $
*/

#ifndef DEVICEMANAGER_CUDA_INCLUDED
#define DEVICEMANAGER_CUDA_INCLUDED

#include <string>

#include "lib/Parameters/commonParameters.h"
#include "lib/IO/bridgeIO.h"

//! Device manger for CUDA.

/*!
   Device manger for CUDA.
   This class is static.
                                    [23 Jul 2023 H.Matsufuru]
 */
class DeviceManager_CUDA {

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
  DeviceManager_CUDA(const DeviceManager_CUDA&);
  DeviceManager_CUDA& operator=(const DeviceManager_CUDA&);

};
#endif
