/*!
        @file    deviceManager_NVSHMEM.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2024-10-29 11:17:40 #$
        @version $LastChangedRevision: 2612 $
*/

#ifndef DEVICEMANAGER_NVSHMEM_INCLUDED
#define DEVICEMANAGER_NVSHMEM_INCLUDED

#include <string>

#include <cuda_runtime.h>

#include "lib/Parameters/commonParameters.h"
#include "lib/IO/bridgeIO.h"

//! Device manger for NVSHMEM.

/*!
   Device manger for NVSHMEM.
   This class is static.
                                    [23 Jul 2023 H.Matsufuru]
 */
class DeviceManager_NVSHMEM {

 public:
  static const std::string class_name;

 private:
  static int m_num_devices;          //!< number of devices.
  static int m_device_id;            //!< device id used in this process.
  static Bridge::VerboseLevel m_vl;  //!< verbose level.

  static int nstream;                //!< number of stream.
  static cudaStream_t* cuda_stream;                                           

 public:
  //! setup: called in main only once.
  static void init(int device_num_offset);

  //! finalization.
  static void finalize();

 private:
  // non-copyable
  DeviceManager_NVSHMEM(const DeviceManager_NVSHMEM&);
  DeviceManager_NVSHMEM& operator=(const DeviceManager_NVSHMEM&);

};
#endif
