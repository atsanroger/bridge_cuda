/*!
      @file    bridgeACC.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2024-07-13 20:19:07 #$
      @version $LastChangedRevision: 2606 $
*/

#ifndef BRIDGEACC_INCLUDED
#define BRIDGEACC_INCLUDED

#include <stdio.h>
#include <cuda_runtime.h>

typedef cudaStream_t cudaStream_t;

//#define MAX_DEVPTR_REF 256
//#define MAX_DEVPTR_REF 1024
// Raised for the batched (MRHS) AMG: solving s=12 propagator columns at once
// keeps s x (block-Krylov basis + V-cycle/smoother/coarse buffers) AFields live
// simultaneously, which exceeds 1024.  The registry deregisters on field
// destruction, so this caps SIMULTANEOUS live fields, not cumulative; 8192 x
// 3 ptrs x 8 B ~ 200 KB, negligible.
#define MAX_DEVPTR_REF 8192

#define CHECK(call)                                    \
  { const cudaError_t error = call;                    \
    if(error != cudaSuccess){                          \
      printf("Error: %s:%d, ", __FILE__, __LINE__);    \
      printf("code:%d, reason: %s\n", error,           \
             cudaGetErrorString(error));               \
      exit(1);                                         \
    }                                                  \
  }

#define CREATE_CUDA_STREAMS(stream_ptr, nstream) \
    do { \
        if ((stream_ptr)  == nullptr) return; \
        *(stream_ptr) = (cudaStream_t*)malloc((nstream) * sizeof(cudaStream_t)); \
        if (*(stream_ptr) == nullptr) return; \
        for (int istream = 0; istream < (nstream); ++istream) { \
            cudaError_t err = cudaStreamCreate(&((*stream_ptr)[istream])); \
            if (err != cudaSuccess) { \
                printf("Failed to create CUDA stream %d: %s\n", istream, cudaGetErrorString(err)); \
                for (int jstream = 0; jstream < istream; ++jstream) { \
                    cudaStreamDestroy((*stream_ptr)[jstream]); \
                } \
                free(*(stream_ptr)); \
                *(stream_ptr) = nullptr; \
                return; \
            } \
        } \
    } while (0)

#define DESTROY_CUDA_STREAMS(stream_ptr, nstream) \
    do { \
        if ((stream_ptr) == nullptr || *(stream_ptr) == nullptr) return; \
        for (int istream = 0; istream < (nstream); ++istream) { \
            cudaStreamDestroy((*stream_ptr)[istream]); \
        } \
        free(*(stream_ptr)); \
        *(stream_ptr) = nullptr; \
    } while (0)

namespace BridgeACC {

void* dev_ptr(void* phost);

void init();
void tidyup();

// real_t = double
void enter_data_create(double* data, size_t size);
void exit_data_delete(double* data, size_t size);

void update_device(double* x, size_t offset, size_t size);
void update_host(double* x, size_t offset, size_t size);

void update_device_asy(double* x, size_t offset, size_t size, cudaStream_t nsteam);
void update_host_asy(double* x, size_t offset, size_t size, cudaStream_t nsteam);

// real_t = float
void enter_data_create(float* data, size_t size);
void exit_data_delete(float* data, size_t size);

void update_device(float* x, size_t offset, size_t size);
void update_host(float* x, size_t offset, size_t size);

void update_device_asy(float* x, size_t offset, size_t size, cudaStream_t nsteam);
void update_host_asy(float* x, size_t offset, size_t size, cudaStream_t nsteam);

}

#endif
