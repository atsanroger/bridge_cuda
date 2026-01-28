/*!
      @file    bridsgeACC.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-05-17 15:25:49 #$
      @version $LastChangedRevision: 2258 $
*/

#include <stdio.h>

#include "inline/define_params.h"
#include "inline/define_index.h"

#include "lib_alt_Accel/BridgeOpenACC/bridgeACC.h"

namespace BridgeACC {

// real_t = double
//====================================================================
void enter_data_create(double* data, const size_t size)
{
#pragma acc enter data create(data[0:size])
}

//====================================================================
void exit_data_delete(double* data, const size_t size)
{
#pragma acc exit data delete(data[0:size])
}

//====================================================================
void update_device(double* data, size_t offset, size_t size)
{
#pragma acc update device(data[offset:size])
}

//====================================================================
void update_host(double* data, size_t offset, size_t size)
{
#pragma acc update host(data[offset:size])
}

// real_t = float
//====================================================================
void enter_data_create(float* data, const size_t size)
{
#pragma acc enter data create(data[0:size])
}

//====================================================================
void exit_data_delete(float* data, const size_t size)
{
#pragma acc exit data delete(data[0:size])
}

//====================================================================
void update_device(float* data, size_t offset, size_t size)
{
#pragma acc update device(data[offset:size])
}

//====================================================================
void update_host(float* data, size_t offset, size_t size)
{
#pragma acc update host(data[offset:size])
}

}
//============================================================END=====          

