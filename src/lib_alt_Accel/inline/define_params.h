/*!
        @file    define_params.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: wlchen $
        @date    $LastChangedDate:: 2025-04-09 14:25:12 #$
        @version $LastChangedRevision: 2535 $
*/

#ifndef ACCEL_DEFINE_PARAMS_INCLUDED
#define ACCEL_DEFINE_PARAMS_INCLUDED

#include <cuda_runtime.h>

#if defined USE_GROUP_SU3
//#include "lib_alt_Accel/inline/define_params_SU3.h"

//temporary setting
#define  NC      3
#define  NVC     6
#define  NDF    18
#define  NDF2    9
#define  ND      4
#define  ND2     2
#define  NCD    12
#define  NVCD   24
#define  NDIM    4

#endif

//#define  NWP     1
#define  NWP     32

#define NUM_WORKERS 1
//#define VECTOR_LENGTH 8
#define VECTOR_LENGTH 64

#define WARP_LENGTH 32
#define MAX_THREAD_PER_BLOCK 1024
#define NS_MAX 64

//#define CEIL_NWP(nst) (nst)
#define CEIL_NWP(nst) (((nst + NWP-1)/NWP) * NWP)
// this definition must be same as ceil_nwp() function below

namespace {

  inline int ceil_nwp(const int nst)
  {
    int size = nst/NWP;
    if( (nst % NWP) > 0) ++size;
    return size * NWP;
  }

  inline int getNumSMs()
  {
    static int numSMs = []()
    {
      cudaDeviceProp prop;
      cudaGetDeviceProperties(&prop, 0);
      return prop.multiProcessorCount;
    }();
    return numSMs;
  }

  
} // namespace

#endif
//============================================================END=====
