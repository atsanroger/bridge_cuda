/*!
        @file    bridge_cuda-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2024-10-29 11:17:40 #$
        @version $LastChangedRevision: 2612 $
*/

#include <stdio.h>

#include <nvshmem.h>
#include <nvshmemx.h>
#include <cuda_runtime.h>

extern void* ptr_host[];
extern void* ptr_dev[];
extern size_t ptr_size[];
extern int devptr_max;

//====================================================================
void enter_data_create(real_t* data, size_t size)
{
  //printf("enter_data_create is called: %d\n", sizeof(real_t));
  //printf("  devptr_max = %d\n", devptr_max);

  if(devptr_max >= MAX_DEVPTR_REF){
    printf("too many device pointer references: %d\n", devptr_max);
    printf("increase MAX_DEVPTR_REF (currently %d)\n", MAX_DEVPTR_REF);
    exit(1);
  }

  ptr_host[devptr_max] = (void*)data;
  //printf("  host pointer = %p\n", ptr_host[devptr_max]);

  size_t count = size * sizeof(real_t);

  ptr_size[devptr_max] = count;

  //printf("  try nvshmem.\n");

  //CHECK(cudaMalloc((void**)&ptr_dev[devptr_max], count));
  void* ptr = nvshmem_malloc(count);

  //printf("  device memory pointr = %p\n", ptr);

  if(ptr != 0){
    ptr_dev[devptr_max] = ptr;
  }else{
    printf("nvshmem_malloc returned null pointer.\n");
    exit(1);
  }
 

  ++devptr_max;

  //printf("enter_data_create finished.\n");

}

//====================================================================
void exit_data_delete(real_t* data, size_t size)
{
  void* dev_mem = dev_ptr((void*)data);
  //real_t* dev_mem = (real_t*)dev_ptr((void*)data);

  //printf("exit_data_delete called.\n");
  //printf("  devptr_max = %d\n", devptr_max);
  //printf("  device memory pointr = %p\n", dev_mem);

  // CHECK(cudaFree((void*)dev_mem));
  nvshmem_free((void*)dev_mem);

  void* ptr = (void*)data;

  int target = -1;
  for(int i = 0; i < devptr_max; ++i){
    if(ptr_host[i] == ptr){
      target = i;
    }
  }

  //printf("  target index = %d\n", target);

  if(target > -1){
    --devptr_max;
    for(int i = target; i < devptr_max; ++i){
      ptr_host[i] = ptr_host[i+1];
      ptr_dev[i]  = ptr_dev[i+1];
      ptr_size[i] = ptr_size[i+1];
    }
    ptr_host[devptr_max] = 0;
    ptr_dev[devptr_max]  = 0;
    ptr_size[devptr_max] = 0;
    return;
  }

  printf("no registered host pointer: %p\n", ptr);
  exit(1);

}

//====================================================================
void update_device(real_t* x, size_t offset, size_t size)
{
  real_t* x_dev = (real_t*)dev_ptr(x);

  //printf("update_device is called: %p\n", x_dev);

  CHECK(cudaMemcpy(&x_dev[offset], &x[offset],
                   size * sizeof(real_t),
                   cudaMemcpyHostToDevice));

}

//====================================================================
void update_host(real_t* x, size_t offset, size_t size)
{
  real_t* x_dev = (real_t*)dev_ptr(x);

  //printf("update_host is called: %p\n", x_dev);

  CHECK(cudaMemcpy(&x[offset], &x_dev[offset],
                   size * sizeof(real_t),
                   cudaMemcpyDeviceToHost));

}
//====================================================================
void update_device_asy(real_t* x, size_t offset, size_t size, cudaStream_t nstream)
{
  real_t* x_dev = (real_t*)dev_ptr(x);

  if (nstream == nullptr) {
      printf("CUDA stream is null, back to synchronous cudaMemcpy.\n");
      update_device(x, offset, size);
      return;
  }

  CHECK(cudaMemcpyAsync(&x_dev[offset], &x[offset],
                   size * sizeof(real_t),
                   cudaMemcpyHostToDevice,
                   nstream));

}

//====================================================================
void update_host_asy(real_t* x, size_t offset, size_t size, cudaStream_t nstream)
{
  real_t* x_dev = (real_t*)dev_ptr(x);

  if (nstream == nullptr) {
      printf("CUDA stream is null, back to synchronous cudaMemcpy.\n");
      update_host(x, offset, size);
      return;
  }

  CHECK(cudaMemcpyAsync(&x[offset], &x_dev[offset],
                   size * sizeof(real_t),
                   cudaMemcpyDeviceToHost,
                   nstream));

}

//====================================================================
