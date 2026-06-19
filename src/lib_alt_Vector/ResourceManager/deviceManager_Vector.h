/*!
        @file    deviceManager_Vector.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
        @version $LastChangedRevision: 2543 $
*/

#ifndef DEVICEMANAGER_VECTOR_INCLUDED
#define DEVICEMANAGER_VECTOR_INCLUDED

#include <string>

#include "lib/Parameters/commonParameters.h"
#include "lib/IO/bridgeIO.h"

//! Device manger for Vector.

/*!
   Device manger for Vector.
   This class is static.
                                    [06 Jun 2019 H.Matsufuru]
 */
class DeviceManager_Vector {

 public:
  static const std::string class_name;

 private:
  static int m_num_devices;          //!< number of devices.
  static int m_device_id;            //!< device id used in this process.
  static Bridge::VerboseLevel m_vl;  //!< verbose level.

 public:
  //! setup: called in main only once.
  static void init();

  //! finalization.
  static void finalize();

 private:
  // non-copyable
  DeviceManager_Vector(const DeviceManager_Vector&);
  DeviceManager_Vector& operator=(const DeviceManager_Vector&);

};
#endif
