/*!
      @file    bridgeACC.cu
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-07-23 18:49:38 #$
      @version $LastChangedRevision: 2534 $
*/

#include "bridgeACC.h"

#include "nvshmem.h"
#include "nvshmemx.h"

namespace BridgeACC {

void*  ptr_host[MAX_DEVPTR_REF];
void*  ptr_dev[MAX_DEVPTR_REF];
size_t ptr_size[MAX_DEVPTR_REF];
int devptr_max = 0;

cudaStream_t cuda_stream[NUM_CUDASTREAM];

//====================================================================
void init()
{
  int nstream = NUM_CUDASTREAM;
  for(int istream = 0; istream < nstream; ++istream){
    cudaError_t err;
    err = cudaStreamCreate(&cuda_stream[istream]);
    if (err != cudaSuccess) {
      printf("Failed to create CUDA stream %d: %s\n",
             istream, cudaGetErrorString(err));
      exit(EXIT_FAILURE);
    }
  }
}

//====================================================================
void tidyup()
{
  int nstream = NUM_CUDASTREAM;
  for (int istream = 0; istream < nstream; ++istream) {
    cudaStreamDestroy(cuda_stream[istream]);
  }

}

//====================================================================
void* dev_ptr(void* phost)
{
  for(int i = 0; i < devptr_max; ++i){
    char* upper = ((char*)ptr_host[i] + ptr_size[i]);
    char* lower = (char*)ptr_host[i];
    if((char*)phost >= lower && (char*)phost < upper){
      size_t offset = (char*)phost - (char*)ptr_host[i];
      return (void*)((char*)ptr_dev[i] + offset);
    }
  }
  printf("no registered device pointer: %p\n", phost);
  exit(1);
  return 0;
}

//====================================================================
__global__
void comm_nvshmem_kernel(char* v2, char* v1, size_t size, int dest)
{
  nvshmem_char_put(v2, v1, size, dest);
}

// API function: char
void comm_nvshmem(char* v2, char* v1, size_t size, int dest, int stream)
{
  char* v2_dev = (char*)dev_ptr(v2);
  char* v1_dev = (char*)dev_ptr(v1);

  comm_nvshmem_kernel<<<1, 1, 0, cuda_stream[stream]>>>(
                                       v2_dev, v1_dev, size, dest);

  nvshmemx_barrier_all_on_stream(cuda_stream[stream]);

}

//====================================================================
void comm_nvshmem(double* v2, double* v1, size_t size, int dest, int stream)
{
  double* v2_dev = (double*)dev_ptr((void*)v2);
  double* v1_dev = (double*)dev_ptr((void*)v1);

  nvshmem_double_put(v2, v1, size, dest);

}

//====================================================================
// API function: float
void comm_nvshmem(float* v2, float* v1, size_t size, int dest, int stream)
{
  char* v2_dev = (char*)dev_ptr((void*)v2);
  char* v1_dev = (char*)dev_ptr((void*)v1);
  size_t size2 = sizeof(float) * size;

  comm_nvshmem_kernel<<<1, 1, 0, cuda_stream[stream]>>>(
                                       v2_dev, v1_dev, size2, dest);

  nvshmemx_barrier_all_on_stream(cuda_stream[stream]);

}


}  // namespace BridgeACC
//============================================================END=====
