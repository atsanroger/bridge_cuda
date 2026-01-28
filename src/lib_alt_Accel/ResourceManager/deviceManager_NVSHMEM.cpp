/*!
        @file    deviceManager_NVSHMEM.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2024-10-29 11:17:40 #$
        @version $LastChangedRevision: 2612 $
*/

#ifdef USE_ACCEL_NVSHMEM

#include "lib_alt_Accel/ResourceManager/deviceManager_NVSHMEM.h"

#include "nvshmem.h"
#include "nvshmemx.h"
#include "mpi.h"

// this branch assumes NVIDIA environment
#include "accel.h"
#include "cuda_runtime.h"

#include "lib/Communicator/communicator.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_Accel/BridgeNVSHMEM/bridgeACC.h"

int DeviceManager_NVSHMEM::m_num_devices = 1;
int DeviceManager_NVSHMEM::m_device_id  = 0;

Bridge::VerboseLevel DeviceManager_NVSHMEM::m_vl = Bridge::CRUCIAL;

const std::string DeviceManager_NVSHMEM::class_name
                                           = "DeviceManager_NVSHMEM";

//====================================================================
void DeviceManager_NVSHMEM::init(int device_num_offset = 0)
{
  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: initialization start.\n",
               class_name.c_str());

  int process_id = Communicator::nodeid();
  vout.general(m_vl, "  process_id = %d\n", process_id);

  vout.general(m_vl, "  NVIDIA environment is used.\n");
  cudaGetDeviceCount(&m_num_devices);
  m_device_id = process_id % m_num_devices + device_num_offset;
  cudaSetDevice(m_device_id);
  vout.general(m_vl, "  number of devices = %d\n", m_num_devices);
  vout.general(m_vl, "  device_id = %d  is used on process = %d\n",
               m_device_id, process_id);

  // nvshmem
  nvshmemx_init_attr_t attr;
  MPI_Comm comm = MPI_COMM_WORLD;
  attr.mpi_comm = &comm;
  nvshmemx_init_attr(NVSHMEMX_INIT_WITH_MPI_COMM, &attr);

  vout.general(m_vl, "  NVSHMEM initialized\n");
  int my_pe = nvshmem_my_pe();
  int n_pes = nvshmem_n_pes();
  vout.general(m_vl, "    my_pe = %d\n", my_pe);
  vout.general(m_vl, "    n_pes = %d\n", n_pes);

  // cudaStream
  BridgeACC::init();
  int nstream = NUM_CUDASTREAM;
  vout.general(m_vl, "  cudaStreams were created%d\n");
  vout.general(m_vl, "  number of streams = %d\n", nstream);

  vout.general(m_vl, "%s: initialization finished.\n",
               class_name.c_str());
}

//====================================================================
void DeviceManager_NVSHMEM::finalize()
{
  vout.paranoiac(m_vl, "%s: finalization.\n", class_name.c_str());

  nvshmem_finalize();

  BridgeACC::tidyup();
}

#endif
//============================================================END=====
