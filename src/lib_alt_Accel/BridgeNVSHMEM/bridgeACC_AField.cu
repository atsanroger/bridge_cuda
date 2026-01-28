/*!
      @file    bridgeACC_AField.cu
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2225 $
*/

#include "inline/define_params.h"
#include "inline/define_index.h"
#include "bridgeACC.h"

#include "bridgeACC_AField.h"


namespace BridgeACC {

//====================================================================
__global__ void copy_kernel(double* v, float* w)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  v[i] = w[i];
}

//====================================================================
void copy(double* y, int nv1, float* x, int nv2, int nin, int nvol)
{
  int size = nin * nvol;

  double* y_dev = (double*)dev_ptr(y);
  float*  x_dev =  (float*)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = size/nth;

  copy_kernel<<<nbl,nth>>>(&y_dev[nv1], &x_dev[nv2]);

  cudaDeviceSynchronize();
}


//====================================================================
__global__ void copy_kernel(float* v, double* w)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  v[i] = w[i];
}

//====================================================================
void copy(float* y, int nv1, double* x, int nv2, int nin, int nvol)
{
  int size = nin * nvol;

  float* y_dev  =  (float*)dev_ptr(y);
  double* x_dev = (double*)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = size/nth;

  copy_kernel<<<nbl,nth>>>(&y_dev[nv1], &x_dev[nv2]);

  cudaDeviceSynchronize();
}

}  // namespace BridgeACC

//============================================================END=====
